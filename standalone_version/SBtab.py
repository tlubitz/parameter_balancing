#!/usr/bin/env python
'''
SBtab
=====
SBtab is a uniforming table format and designed for the use in
Systems Biology. Furthermore, it is a useful format to import stored
information into Python objects to process and manipulate it.

See specification for further information.
'''
import re
import copy
import tablib
try:
    from . import tablibIO
    from . import misc
except:
    import tablibIO
    import misc


class SBtabError(Exception):
    '''
    Base class for errors in the SBtab class.
    '''
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message


class SBtabTable():
    '''
    SBtabTable (version 0.9.0 06/10/2010)
    '''
    def __init__(self, table_string, filename):
        '''
        Creates SBtab Python object from tablib object.

        Parameters
        ----------
        table : str
            Containing one SBtab, uncut. Directly from f.read().
        filename : str
            Filename with extension.
        '''
        # Needed to be able to adress it from outside of the class for
        # writing and reading
        self.filename = filename
        self.delimiter = misc.check_delimiter(table_string)
        table_list = self.cut_table_string(table_string)

        # check if ascii stuff is violated
        try: (self.table, self.str_tab) = self.check_ascii(table_list)
        except: raise SBtabError('''This is not a valid SBtab file. Try to
        check your file with the SBtab validator or read the
        SBtab specification.''')

        # Delete tablib header to avoid complications
        if self.table.headers: self.table.headers = None

        # Create all necessary variables
        self.initialize_table()
        self.sbtab_dataset = []
        self.tables_without_name = []

    def cut_table_string(self, table_string):
        '''
        the SBtab is initially given as one long string;
        cut down string into list to harvest content
        '''
        table_list = []
        for row in table_string.split('\n'):
            table_list.append(row.split(self.delimiter))
        return table_list

    def return_table_string(self):
        '''
        sometimes the file is required as a string (e. g. for
        writing files to harddisk; return string
        '''
        table_string = [self.header_row]
        table_string.append('\t'.join(self.columns))
        for row in self.value_rows:
            table_string.append('\t'.join(row))
        return '\n'.join(table_string)

    def initialize_table(self):
        '''
        Loads table informations and class variables.
        '''
        # Read the header row from table
        self.header_row = self._get_header_row()

        # Read the table information from header row
        (self.table_type,
         self.table_name,
         self.table_document,
         self.table_version) = self.get_table_information()

        # Read the columns of the table
        (self.columns, self.columns_dict) = self.get_columns()

        self.delimiter = misc.check_delimiter('\n'.join(self.str_tab))

        # Read data rows
        self.value_rows = self.get_rows(self.table_type)

        # Update the list and tablib object
        self.update()

    def check_ascii(self, table):
        '''
        Checks for ASCII violations, so that the parser will not crash
        if these occur.

        Parameters
        ----------
        table : list of str
            SBtab table as a list of row strings.
        '''
        new_table = []

        for row in table:
            new_row = []
            for entry in row:
                try:
                    new_row.append(str(entry).strip())
                except:
                    new_row.append('''Ascii violation error!
                                      Please check input file!''')
            new_table.append('\t'.join(new_row))

        tablibtable = tablibIO.importSetNew('\n'.join(new_table),
                                            self.filename + '.csv')

        return tablibtable, new_table

    def _get_header_row(self):
        '''
        Extracts the declaration row from the SBtab file.
        '''
        header_row = None
        # Find header row
        for row in self.table:
            for entry in row:
                if str(entry).startswith('!!'):
                    header_row = row
                    break
                elif str(entry).startswith('"!!'):
                    rm1 = row.replace('""', '#')
                    rm2 = row.remove('"')
                    header_row = rm2.replace('#', '"')
                    break

        # Save string or raise error
        if not header_row:
            raise SBtabError('''This is not a valid SBtab table, please use
            validator to check format or have a look in the specification!''')
        else:
            header_row = ' '.join(header_row)

        # Replace double quotes by single quotes
        stupid_quotes = ['"', '\xe2\x80\x9d', '\xe2\x80\x98', '\xe2\x80\x99',
                         '\xe2\x80\x9b', '\xe2\x80\x9c', '\xe2\x80\x9f',
                         '\xe2\x80\xb2', '\xe2\x80\xb3', '\xe2\x80\xb4',
                         '\xe2\x80\xb5', '\xe2\x80\xb6', '\xe2\x80\xb7']

        for squote in stupid_quotes:
            try: header_row = header_row.replace(squote, "'")
            except: pass

        # Split header row
        header_row = header_row.split(' ')

        # Delete spaces in header row
        while '' in header_row:
            header_row.remove('')

        header = ""
        for x in header_row[:-1]:
            header += x + ' '
        header += header_row[-1]

        return header

    def get_table_information(self):
        '''
        Reads declaration row and stores the SBtab table attributes.
        '''
        no_name_counter = 0

        # Save table type, otherwise raise error
        try: table_type = self.get_custom_table_information('TableType')
        except: raise SBtabError('The TableType of the SBtab is not defined!')

        # Save table name, otherwise give name with number of unnamed tables
        try: table_name = self.get_custom_table_information('TableName')
        except:
            self.tables_without_name.append(table_type)
            for table_no_name in self.tables_without_name:
                if table_type == table_no_name:
                    no_name_counter = no_name_counter + 1
            table_name = table_type.capitalize() + '_' + str(no_name_counter)
            self.header_row += " TableName='" + table_name + "'"

        # Save table document, otherwise return None
        try: table_document = self.get_custom_table_information('Document')
        except: table_document = None

        # save table version, otherwise return None
        try: table_version = self.get_custom_table_information('SBtabVersion')
        except: table_version = None

        return table_type, table_name, table_document, table_version

    def get_custom_table_information(self, attribute_name):
        '''
        Retrieves the value of a table attribute in the declaration line

        Parameters
        ----------
        attribute_name : str
           Name of the table attribute.
        '''
        if re.search("%s='([^']*)'" % attribute_name,
                     self.header_row) is not None:
            return re.search("%s='([^']*)'" % attribute_name,
                             self.header_row).group(1)
        else:
            raise SBtabError('''The %s of the SBtab is
                                not defined!''' % attribute_name)

    def get_columns(self):
        '''
        Extract column headers, add mandatory first column name if necessary.
        '''
        # Save list of main columns
        for row in self.table:
            for entry in row:
                if str(row[0]).startswith('!') \
                   and not str(row[0]).startswith('!!'):
                    column_names = list(row)
                    break

        # Get column positions
        columns = {}
        for i, column in enumerate(column_names):
            columns[column] = i

        return column_names, columns

    def get_rows(self, table_type='table', inserted=False):
        '''
        Extract the rows of the SBtab, add first column if necessary.

        Parameters
        ----------
        table_type : str (default 'table')
            Attribute TableType of the SBtab table. Only necessary, if first
            column was set automatically.
        inserted : Boolean
            True, if mandatory first column was set automatically.
        '''
        # Add row to list value_rows
        # if row doesn't contain entries starting with '!'
        value_rows = []

        # Add to comments, if row starts with '%'
        self.comments = []

        for row in self.table:
            if str(row[0]).startswith('!'):
                continue
            for i, entry in enumerate(row):
                if str(entry).startswith('%'):
                    self.comments.append(list(row))
                    break
                else:
                    if len(row) == i + 1:
                        value_rows.append(list(row))

        # Insert value column if mandatory column was added automatically
        if inserted:
            if table_type == 'table':
                for i, row in enumerate(value_rows):
                    row.insert(0, 'TableRow_' + str(i + 1))
            else:
                for i, row in enumerate(value_rows):
                    row.insert(0, table_type[0].upper() +
                               table_type[-1].lower() + str(i + 1))

        return value_rows

    def change_value(self, row, column, new):
        '''
        Change single value in the SBtab table by position in the table.

        Parameters
        ----------
        row : int
            Number of rows in the table. First row is number 1.
        column : int
            Number of columns in the table. First column is number 1.
        new : str
            New entry.
        '''
        self.value_rows[row - 1][column - 1] = new

        # Update object
        self.update()

    def change_value_by_name(self, name, column_name, new):
        '''
        Change singe value in the SBtab by name of column
        and of the first row entry.

        Parameters
        ----------
        row : str
            Name of the entry in the first column.
        column : str
            Name of the column (without '!')
        new : str
            New entry.
        '''
        col = self.columns_dict['!' + column]
        for r in self.value_rows:
            if r[0] == name:
                r[col] = new

        # Update object
        self.update()

    def create_list(self):
        '''
        Creates a list object of the SBtab Python object.
        '''
        # Create new list
        sbtab_list = []

        # Append the parts header row, main column row and value
        # rows to the list
        sbtab_list.append(self.header_row)
        sbtab_list.append(self.columns)
        sbtab_list.append(self.value_rows)

        return sbtab_list

    def create_dataset(self):
        '''
        Creates a tablib object of the SBtab Python object.
        '''
        # Initialize empty variables for conversion
        sbtab_temp = []
        self.sbtab_dataset = tablib.Dataset()

        # Create list of header
        header = [self.header_row]

        # Delete spaces in header, main column and data rows
        header = [x.strip(' ') for x in header]
        self.columns = [x.strip(' ') for x in self.columns]
        for row in self.value_rows:
            try:
                for entry in row:
                    entry = entry.strip(' ')
            except:
                continue

        # Add header, main column and data rows to temporary list object
        sbtab_temp.append(header)
        sbtab_temp.append(self.columns)
        for row in self.value_rows:
            sbtab_temp.append(row)

        # Delete all empty entries at the end of the rows
        sb1 = []
        for row in sbtab_temp:
            if row[0] != '':
                sb1.append(row)

        # Make all rows the same length
        longest = max([len(x) for x in sb1])

        for row in sb1:
            if len(row) < longest:
                for i in range(longest - len(row)):
                    row.append('')
                self.sbtab_dataset.append(row)
            else:
                self.sbtab_dataset.append(row)

        return self.sbtab_dataset

    def add_row(self, row_list, position=None):
        '''
        Adds row to the table, if postion is None at the end of it.

        Parameters
        ----------
        row_list : list
            List of strings, containing the entries of the new row.
        position : int
            Position of new row in the table, 0 is on top.
        '''
        # Empty column to fill up sbtab_dataset with ''
        empty_list = []

        # Create temporary work copy
        sbtab_dataset = self.table

        # If new row is to small, add empty entries to new row
        if len(row_list) < len(sbtab_dataset.dict[0]):
            for i in range(len(sbtab_dataset.dict[0]) - len(row_list)):
                row_list.append('')

        # If new row is too long, add empty entries to sbtab_dataset
        elif len(row_list) > len(sbtab_dataset.dict[0]):
            for i in range(len(sbtab_dataset.dict[0])):
                empty_list.append('')
            for i in range(len(row_list) - len(sbtab_dataset.dict[0])):
                sbtab_dataset.rpush_col(empty_list)

        # If no position is set, add new row to the end
        if position is None:
            sbtab_dataset.rpush(row_list)
        else:
            sbtab_dataset.insert(position, row_list)

        # Update object
        self.table = sbtab_dataset
        self.initialize_table()

    def remove_row(self, position):
        '''
        Removes row from the table

        Parameters
        ----------
        position : int
            Position of row to be removed. Starting with 1.
        '''
        # Create temporary work copy
        sbtab_dataset = self.table

        del sbtab_dataset[position + 1]

        # Update object
        self.table = sbtab_dataset
        self.initialize_table()

    def add_column(self, column_list, position=None):
        '''
        Adds a column to the table, if position is None at the end of it.

        Parameters
        ----------
        column_list : list
            List of strings, containing the entries of the new column.
        position : int
            Position of new column in the table, 0 is right.
        '''
        # Empty column to fill up sbtab_dataset with ''
        empty_list = []

        # If new column is to small, add empty entries to new column
        if len(column_list) < (len(self.sbtab_dataset.dict) - 1):
            for i in range((len(self.sbtab_dataset.dict) - 1) -
                           len(column_list)):
                column_list.append('')

        # If new column is to long, add empty entries to sbtab_dataset
        elif len(column_list) > (len(self.sbtab_dataset.dict) - 1):
            for i in range(len(self.sbtab_dataset.dict[0])):
                empty_list.append('')
            for i in range(len(column_list) -
                           (len(self.sbtab_dataset.dict) - 1)):
                self.value_rows.append(empty_list)
                empty_list = copy.deepcopy(empty_list)

        # If no position is set, add new column to the end
        if not position:
            for i, row in enumerate(self.value_rows):
                row.append(column_list[i + 1])
            self.columns_dict[column_list[0]] = len(self.columns)
            self.columns = self.columns_dict.keys()
        else:
            for i, row in enumerate(self.value_rows):
                row.insert(position - 1, column_list[i + 1])
            self.columns_dict[column_list[0]] = position - 1
            self.columns = self.columns_dict.keys()

        # Update object
        self.update()

    def remove_column(self, position):
        '''
        Removes column from the table.

        Parameters
        ----------
        position : int
            Position of column to be removed. Sarting with 1.
        '''
        # Remove entries on position
        for row in self.value_rows:
            del row[position + 1]
        for column in self.columns_dict.keys():
            if self.columns_dict[column] == position - 1:
                del self.columns_dict[column]

        # Update object
        self.update()

    def write_sbtab(self, format_type, filename=None):
        '''
        Writes SBtab tablib object to file.

        Parameters
        ----------
        format_type : str
            File extension of the SBtab file. ('tsv', 'csv', 'tab', 'xls')
        filename : str
            Filename of the SBtab file without extension. Default is filename.
        '''
        if not filename:
            filename = self.filename[:-4]
        if format_type == 'tsv' or format_type == 'tab':
            tablibIO.writeTSV(self.sbtab_dataset, filename)
        elif format_type == 'csv':
            tablibIO.writeCSV(self.sbtab_dataset, filename)
        elif format_type == 'ods':
            tablibIO.writeODS(self.sbtab_dataset, filename)
        elif format_type == 'xls':
            tablibIO.writeXLS(self.sbtab_dataset, filename)
        else:
            raise SBtabError('''The given file format is not supported: %s.
                                Please use ".tsv", ".csv", ".tab"
                                or ".xls" instead.''' % (format_type))

    def duplicate(self):
        '''
        Creates a copy of the SBtab object.
        '''
        sbtab = copy.deepcopy(self)

        return sbtab

    def update(self):
        '''
        Updates the SBtab instance, list object, and tablib dataset.
        '''
        # Create tablib Dataset instance with new SBtab table
        self.table = self.create_dataset()
        # Create list instance with new SBtab table
        self.sbtab_list = self.create_list()

    def create_sbtab_dict(self):
        '''
        Creates a dict instance of the SBtab table.
        Keys are the column names, values are dicts. These contain the entries
        of the table. Keys are the entries in the first column, values are the
        current entries in the certain column.
        '''
        sbtab_dicts = {}
        for column_name in self.columns:
            sbtab_dicts[column_name] = {}
            for row in self.value_rows:
                sbtab_dicts[column_name][row[0]] = row[self.columns_dict[column_name]]

        return sbtab_dicts

    def transpose_table(self):
        '''
        Transposes SBtab table. Switches columns and rows.
        '''
        # Initialize new table data
        trans_columns = []
        trans_columns_dict = {}
        trans_value_rows = []

        # Save old table data
        columns = self.columns
        value_rows = self.value_rows

        # Append first entry to new column
        trans_columns.append(columns.pop(0))

        # Set new rows
        for column in columns:
                trans_value_rows.append([column])

        # Set new values in tables
        for row in value_rows:
            trans_columns.append(row.pop(0))
            for i, entry in enumerate(row):
                trans_value_rows[i].append(entry)

        # Write new columns dict
        for i, column in enumerate(trans_columns):
            trans_columns_dict[column] = i

        # Overwrite old table data
        self.columns = trans_columns
        self.columns_dict = trans_columns_dict
        self.value_rows = trans_value_rows

        self.update()

    def to_data_frame(self):
        import pandas as pd
        column_names = map(lambda s: s[1:], self.columns)
        df = pd.DataFrame(data=self.get_rows(), columns=column_names)
        return df
