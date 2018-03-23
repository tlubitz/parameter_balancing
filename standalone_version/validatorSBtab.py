"""
SBtab Validator
===============

Python script that validates SBtab files

See specification for further information.
"""
#!/usr/bin/env python
from . import SBtab
from . import tablibIO
import re
import collections
import sys
from . import misc

class SBtabError(Exception):
    '''
    Base class for errors in the SBtab validation class.
    '''
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message


class ValidateTable:
    '''
    Validator (version 0.9 06/10/2015).
    Checks SBtab file and SBtab object.
    '''
    def __init__(self, table, sbtab_name, def_table=None,def_name=None):
        '''
        Initialises validator and starts check for file and table format.

        Parameters
        ----------
        table : str
            SBtab file as a string.
        sbtab_name : str
            File path of the SBtab file.
        def_table : str
            SBtab definition table as string representation.
        def_name : str
            SBtab definition table name.
        '''
        delimiter    = misc.check_delimiter(table)
        sbtab_tablib = tablibIO.importSetNew(table,sbtab_name,delimiter)
        
        if not def_name:
            def_name = 'definitions.tsv'
        if not def_table:
            try:
                default_def = open(def_name,'r')
                def_table   = default_def.read()
                default_def.close()
            except:
                print('Definition file could not be loaded, so the validation could not be started. Please provide definition file as argument or make it is located in the same directory as this script.')
                sys.exit()
        
        # import definitions from definition table
        definition_table = tablibIO.importSetNew(def_table,def_name,separator='\t')
        definition_sbtab = SBtab.SBtabTable(definition_table, def_name)
        self.definitions = definition_sbtab.sbtab_list

        # create set of valid table types
        self.allowed_table_types = list(set([row[2] for row in self.definitions[2:][0]]))

        # create dict of valid column names per table type
        self.allowed_columns = {}
        for table_type in self.allowed_table_types:
            self.allowed_columns[table_type] = [row[0] for row in self.definitions[2:][0] if row[2] == table_type]

        # initialize warning string
        self.warnings = []
        # define self variables
        self.table = sbtab_tablib
        self.filename = sbtab_name

        # check file format and header row
        self.checkTableFormat()

        # try creating SBtab instance
      
        self.sbtab = SBtab.SBtabTable(self.table, self.filename)

        self.column2format = {}
        defs = self.definitions[2]
        for row in defs:
            if row[3] == self.sbtab.table_type:
                self.column2format[row[1]] = row[4]

        # remove empty column headers
        f_columns = []
        for element in self.sbtab.columns:
            if element == '': pass
            else: f_columns.append(element)
        self.sbtab.columns = f_columns

        # determine headers
        self.determineHeaders()

        # check SBtab object for validity
        self.checkTable()

        '''
        if len(self.warnings) > 0:
            print self.warnings
        else:
            print 'No warnings detected!'
        '''
        
    def checkTableFormat(self):
        '''
        Validates format of SBtab file, checks file format and header row.
        '''
        # Check tablib heade

        if self.table.headers:
            self.warnings.append('Warning: Tablib header is set, will be removed. This feature is not supported.')
            self.table.headers = None
        # save table rows in variable
        self.rows_file = self.table.dict

        # save header row
        header_row = self.rows_file[0]
        while '' in header_row:
            header_row.remove('')
        header = ""
        for x in header_row[:-1]:
            header += x + ' '
        header += header_row[-1]

        #Replace double quotes by single quotes
        stupid_quotes = ['"','\xe2\x80\x9d','\xe2\x80\x98','\xe2\x80\x99','\xe2\x80\x9b','\xe2\x80\x9c','\xe2\x80\x9f','\xe2\x80\xb2','\xe2\x80\xb3','\xe2\x80\xb4','\xe2\x80\xb5','\xe2\x80\xb6','\xe2\x80\xb7']

        for squote in stupid_quotes:
            try: header = header.replace(squote, "'")
            except: pass


        # check for valid header row
        if not header.startswith('!!'):
            self.warnings.append('Error: The header row of the table does not start with "!!SBtab". This file cannot be validated.')
        if not re.search("TableType='([^']*)'", header):
            self.warnings.append('Error: The attribute TableType is not defined in the SBtab table; This file cannot be validated.')
        if not re.search("TableName='([^']*)'", header):
            self.warnings.append('Warning: The (optional) attribute TableName is not defined in the SBtab table.')

        # check for possible table content and main columns
        columns_row = self.rows_file[1]
        columns = ""
        for x in columns_row[:-1]:
            columns += x + ' '
        columns += columns_row[-1]
        # check length of table
        if len(self.rows_file) < 3:
            self.warnings.append('Warning: The table contains no information: ' + \
                header)
        else:
            self.main_column_count = None
            # checks for existing main column
            if not columns.startswith('!'):
                self.warnings.append('Warning: The main column row of the table does not start with "!": ' + \
                    columns)
               
    def determineHeaders(self):
        '''
        determines the headers and saves their position
        '''
        self.id_column=None
        for i,column in enumerate(self.sbtab.columns):
            if column == '!ID':
                self.id_column = i

    def checkTable(self):
        '''
        Validates the mandatory format of the SBtab in accordance to the TableType attribute.
        '''
        column_check = True
        # general stuff
        # 1st: check validity of table_type and save table type for later tests
        if not self.sbtab.table_type in self.allowed_table_types:
            self.warnings.append('Warning: The SBtab file has an invalid TableType in its header: ' + self.sbtab.table_type)
            column_check = False

        # 2nd: check the important first column
        #if self.sbtab.unique_key == 'True':
        #    first_column = '!' + self.sbtab.table_type
        #    if not self.rows_file[1][0] == first_column:
        #        self.warnings.append('The first column of the file does not correspond with the given TableType ' + \
        #                             self.sbtab.table_type + ' and will be filled automatically.')

        #New 2nd: check for ID column
        #if not '!ID' in self.sbtab.columns_dict:
        #    warning = 'There is no !ID column in your SBtab file. This may compromise automated parsing strongly. Please add this column.'
        #    self.warnings.append(warning) 

        unique = []
        # 2,5nd: very important: check if the identifiers start with a digit; this is not allowed in SBML!
        for row in self.sbtab.value_rows:
            if not self.id_column: break
            identifier = row[self.id_column]
            if not identifier in unique: unique.append(identifier)
            else:
                warning = 'Warning: There is an identifier that is not unique. Please change that: '+str(identifier)
                self.warnings.append(warning)
            try:
                int(identifier[0])
                self.warnings.append('Warning: There is an identifier that starts with a digit; this is not permitted for the SBML conversion: '+identifier)
            except: pass

        #not required for our parameter balancing purposes; rather misleading
        #if not '!Name' in self.sbtab.columns_dict:
        #    ident = False
        #    for it in self.sbtab.columns_dict.keys():
        #        if it.startswith('!Identifier'): ident = True
        #    if not ident:
        #        warning = 'Please use at least a !Name or an !Identifier column to better characterise the entries of your SBtab file.'
        #        self.warnings.append(warning)

                
        # Check if there is at least a SumFormula or an identifier to characterise a Reaction
        if self.sbtab.table_type == 'Reaction':
            if not '!ReactionFormula' in self.sbtab.columns_dict.keys():
                for it in self.sbtab.columns_dict.keys():
                    ident = False
                    if it.startswith('!Identifier'): ident = True
                if not ident:
                    warning = 'Error: A Reaction SBtab needs at least a column !ReactionFormula or an !Identifier column to be characterised.'
                    self.warnings.append(warning)

        if self.sbtab.table_type == 'Quantity':
            if not '!Unit' in self.sbtab.columns_dict.keys():
                warning = 'Error: A Quantity SBtab requires the column "'"Unit"'". Please add this column to the SBtab file.'
                self.warnings.append(warning)

        # 3rd: check the validity of the given column names
        if column_check:
            for column in self.sbtab.columns:
                if not column.replace('!', '') in self.allowed_columns[self.sbtab.table_type] and not ('Identifiers:') in column and not ('ID:urn.') in column:
                    self.warnings.append('Warning: The SBtab file has an unknown column: ' + \
                        column + '.\n \t Please use only supported column types!')
            #if not self.sbtab.columns_dict['!' + self.sbtab.table_type] == 0:
            #    self.warnings.append('The SBtab primary column is at a wrong position.')
        else:
            self.warnings.append('Error: The SBtab TableType is "unknown", therefore the main columns cannot be checked!')

        # 4th: check the length of the different rows
        for row in self.sbtab.value_rows:
            # check the content of the main column (first one) for empty entries
            #if row[self.sbtab.columns_dict['!' + self.sbtab.table_type]] == '':
            #    self.warnings.append('The SBtab includes a row with an undefined identifier in the row: \n' + str(row))
            # check the rows for entries starting with + or -
            if str(row[0]).startswith('+') or str(row[0]).startswith('-'):
                self.warnings.append('Warning: An identifier for a data row must not begin with "+" or "-": \n' + str(row))
            if '!ReactionFormula' in self.sbtab.columns_dict:
                if not '<=>' in row[self.sbtab.columns_dict['!ReactionFormula']]:
                    warning = 'There is a sum formula that does not adhere to the sum formula syntax from the SBtab specification: '+str(row[self.sbtab.columns_dict['!ReactionFormula']])
                    self.warnings.append(warning)
                # check the rows for entries containing . or :     OMITTED; should be solved, if the user puts the entryfield in quotes "
                #if ',' in list(row[self.sbtab.columns_dict[column]]):
                #    self.warnings.append('A data row must not include commas, but this one does: \n' + \
                #        str(row[self.sbtab.columns_dict[column]]))
                # raise SBtabError('An identifier for a data row must not
                # include ":" or ".": \n'+str(row))
            for i,it in enumerate(row):
                if it == '': continue
                if self.sbtab.columns[i][1:].startswith('Identifier'): req_format = 'string'
                else:
                    try: req_format = self.column2format[self.sbtab.columns[i][1:]]
                    except: continue
                if req_format == 'Boolean':
                    if it != 'True' and it != 'False' and it != 'TRUE' and it != 'FALSE' and it != '0' and it != '1':
                        warning = 'Warning: The column '+self.sbtab.columns[i][1:]+' holds a value that does not conform with the assigned column format '+req_format+': '+it
                        self.warnings.append(warning)
                elif req_format == 'float':
                    try: float(it)
                    except:
                        warning = 'Warning: The column '+self.sbtab.columns[i][1:]+' holds a value that does not conform with the assigned column format '+req_format+': '+it
                        self.warnings.append(warning)
                elif req_format == '{+,-,0}':
                    if it != '+' and it != '-' and it !='0':
                        warning = 'Warning: The column '+self.sbtab.columns[i][1:]+' holds a value that does not conform with the assigned column format '+req_format+': '+it
                        self.warnings.append(warning)

        # 5th: are there duplicate columns?
        for column in collections.Counter(self.sbtab.columns).items():
            if column[1] > 1:
                self.warnings.append('Warning: There was a duplicate column in this SBtab file. Please remove it: ' + str(column[0]))
        

    def returnOutput(self):
        '''
        Returns the warnings from the validation process.
        '''
        return self.warnings


class ValidateFile:
    '''
    Validates file and checks for valid format.
    '''
    def __init__(self, sbtab_file, filename):
        '''
        Initialisation of file validation.

        Parameters
        ----------
        sbtab_file : str
           SBtab file as string representation.
        filename : str
           SBtab file name.
        '''        
        # initialize warning string
        self.warnings = []

        self.file = sbtab_file
        self.filename = filename

        self.validateExtension()
        #self.validateFile(self.file)

    def validateExtension(self):
        '''
        Checks the extension of the file for invalid formats.
        '''
        valid_extensions = ['tsv','csv','xls','tab']
        if not self.filename[-3:] in valid_extensions:
            self.warnings.append('Error: The file extension is not valid for an SBtab file.')


    def validateFile(self, sbtab_file):
        '''
        Validates file format and checks for possible problems.

        Parameters
        ----------
        sbtab_file : str
           SBtab file as string representation.
        '''
        rows = []
        for line in sbtab_file:
            rows.append(line)

        length = len(rows[0])
        for i, row in enumerate(rows):
            if not row:
                self.warnings.append('Warning: The file contains an empty row in line: ' + str(i))
            if not len(row) == length:
                self.warnings.append('Warning: The lengths of the rows are not identical.\n This will be adjusted automatically.')
            length = len(row)

    def checkSeparator(self):
        '''
        Finds the separator of the file.
        '''
        sep = False

        for row in self.file.split('\n'):
            if row.startswith('!!'): continue
            if row.startswith('!'):
                s = re.search('(.)(!)',row[1:])
                #if there is only one column, we have to define a default separator.
                #let's use a tab.
                try: sep = s.group(1)
                except: sep = '\t'

        return sep

    def returnOutput(self):
        '''
        Returns the warnings from the validation process.
        '''
        return self.warnings

if __name__ == '__main__':

    try: sys.argv[1]
    except:
        print('You have not provided input arguments. Please start the script by also providing an SBtab file and the required definition file: >python validatorSBtab.py SBtab.csv definition.tsv')
        sys.exit()

    file_name    = sys.argv[1]
    sbtab_file_o = open(file_name,'r')
    sbtab_file   = sbtab_file_o.read()
    sbtab_file_o.close()
    
    delimiter    = misc.getDelimiter(sbtab_file)
    sbtab_tablib = tablibIO.importSetNew(sbtab_file,file_name,delimiter)

    try:
        default_def = sys.argv[2]
        def_file    = open(default_def,'r')
        def_tab = def_file.read()
        def_file.close()
    except:
        def_tab = None

    validator_output     = []
    Validate_file_class  = ValidateFile(sbtab_file,file_name)
    validator_output.append(Validate_file_class.returnOutput())
    Validate_table_class = ValidateTable(sbtab_file,file_name,def_tab)
    #Validate_table_class = ValidateTable(sbtab_tablib,file_name,def_tab)
    validator_output.append(Validate_table_class.returnOutput())

    warned = False
    for warning in validator_output:
        if warning != []:
            print('WARNINGS: ',warning)
            warned = True

    if not warned:
        print('The SBtab file is valid.')
            



