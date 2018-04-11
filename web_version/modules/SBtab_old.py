#!/usr/bin/env python
import re

class SBtabError(Exception):
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return self.message

class SBtabTable():
    '''
    SBtab Table
    '''
    def __init__(self,table,filename,table_type=None):
        '''
        initialize the SBtab table
        '''
        self.table_rows = []
        for row in table:
            if not row == ['']:
                self.table_rows.append(row)

        self.separator = '\t'

        if table_type:
            self.table_type = table_type
        else: self.getHeaderRow()
        
        self.getColumns()
        self.getColumnProperties()
        try: self.getColumnProperties()
        except: raise SBtabError('The specification row of the SBtab is invalid (see example files again).')
        self.getRows()
        self.initializeColumns()

    def makeExMarks(self):
        new_rows = []
        for row in self.table_rows:
            if row.startswith('QuantityType'):
                old_column_row = row.split(self.separator)
                new_column_row = []
                for item in old_column_row:
                    new_column_row.append('!'+item)
                new_rows.append(self.separator.join(new_column_row))
            else:
                new_rows.append(row)
        self.table_rows = new_rows
        
    def getHeaderRow(self):
        '''
        extracts the !!-header row from the SBtab file and its information
        '''
        self.header_row = None
       
        for row in self.table_rows:
            if row.startswith('!!'):
                self.header_row = row

        try: self.table_type = re.search('TableType="([^"]*)"',self.header_row).group(1)
        except: self.table_type = None        #will raise Error in validator

        try: self.table_document = re.search('Document="([^"]*)"',self.header_row).group(1)
        except: self.sbtab_document = None

        try: self.table_level = re.search('Level="([^"]*)"',self.header_row).group(1)
        except: self.sbtab_level = None
        
        try: self.table_version = re.search('Version="([^"]*)"',self.header_row).group(1)
        except: self.sbtab_version = None

        try: self.table_name = re.search('Table="([^"]*)"',self.header_row).group(1)
        except: self.sbtab_name = None

    def getColumns(self):
        '''
        extract the column names of the Reaction SBtab
        '''
        #main_name = '!'+self.table_type.capitalize()

        for row in self.table_rows:
            if row.startswith('!') and not row.startswith('!!'):
                self.column_names = row.split(self.separator)
                break
            elif row.startswith('QuantityType'):
                self.column_names = row.split(self.separator)
                break

    def getColumnProperties(self):
        '''
        extract the subcolumns of the Reaction SBtab
        '''
        main_name                 = '!'+self.table_type.capitalize()
        self.column_property_rows = []
        for row in self.table_rows:
            if row.startswith('!') and not row.startswith(main_name) and not row.startswith('!!'):
                self.column_property_rows.append(row.split(self.separator))
            else: break

    def getRows(self):
        '''
        extract the rows of the Reaction SBtab
        '''
        self.value_rows = []
        for row in self.table_rows:
            split_row = row.split(self.separator)
            if not row.startswith('!') and not row.startswith(' ') and not row == '':
                #if len(split_row) == len(self.column_names):
                self.value_rows.append(split_row)

    def initializeColumns(self):
        '''
        initializes the column indices
        '''
        typeFunction = 'self.initializeColumns'+self.table_type+'()'
        eval(typeFunction)
        #try: eval(typeFunction)
        #except: raise SBtabError('The SBtab TableType is invalid: '+self.table_type)

    def initializeColumnsReaction(self):
        '''
        initialize specific columns for the SBtab type Reaction
        '''
        self.main_column       = None
        self.r_column          = None
        self.name_column       = None
        self.sbmlid_column     = None
        self.sumformula_column = None
        self.loc_column        = None
        self.enz_column        = None
        self.kinlaw_column     = None
        self.enzs_column       = None
        self.enzp_column       = None
        self.metreg_column     = None
        self.br_column         = None
        self.be_column         = None
        self.miriam_column     = None

        for i,column in enumerate(self.column_names):
            if column == '!Reaction':
                self.r_column    = i
                self.main_column = i
            elif column == '!Name': self.name_column = i
            elif column == '!SBML::reaction::id': self.sbmlid_column = i
            elif column == '!SumFormula': self.sumformula_column = i
            elif column == '!Location': self.loc_column = i
            elif column == '!Enzyme': self.enz_column = i
            elif column == '!KineticLaw': self.kinlaw_column = i
            elif column == '!Enzyme SBML::species::id': self.enzs_column = i
            elif column == '!Enzyme SBML::parameter::id': self.enzp_column = i
            elif column == '!MetabolicRegulators': self.metreg_column = i
            elif column == '!BuildReaction': self.br_column = i
            elif column == '!BuildEnzyme': self.be_column = i
            elif column.startswith('!MiriamID'): self.miriam_column = i

    def initializeColumnsCompound(self):
        '''
        initialize specific columns for the SBtab type Compound
        '''
        self.main_column   = None
        self.c_column      = None
        self.name_column   = None
        self.sbmlid_column = None
        self.sptype_column = None
        self.loc_column    = None
        self.charge_column = None
        self.const_column  = None
        self.miriam_column = None

        for i,column in enumerate(self.column_names):
            if column == '!Compound':
                self.main_column = i
                self.c_column    = i
            elif column == '!Name': self.name_column = i
            elif column == '!SBML::species::id': self.sbmlid_column = i
            elif column == '!SBML::speciestype::id': self.sptype_column = i
            elif column == '!SpeciesType': self.sptype_column = i
            elif column == '!Location': self.loc_column = i
            elif column == '!Charge': self.charge_column = i
            elif column == '!Constant': self.const_column = i
            elif column.startswith('!MiriamID'): self.miriam_column = i

    def initializeColumnsEnzyme(self):
        '''
        initialize specific columns for the SBtab type Enzyme
        '''
        self.main_column   = None
        self.e_column      = None
        self.name_column   = None
        self.r_column      = None
        self.kinlaw_column = None
        self.metreg_column = None    #needed?
        self.gene_column   = None    #needed?
        self.genes_column  = None    #needed?
        self.miriam_column = None

        for i,column in enumerate(self.column_names):
            if column == '!Enzyme':
                self.main_column = i
                self.e_column    = i
            elif column == '!Name': self.name_column = i
            elif column == '!CatalysedReaction': self.r_column = i
            elif column == '!KineticLaw': self.kinlaw_column = i
            elif column == '!MetabolicRegulators': self.metreg_column = i
            elif column == '!Gene': self.gene_column = i
            elif column == '!GeneBooleanFormula': self.genes_column = i
            elif column.startswith('!MiriamID'): self.miriam_column = i

    def initializeColumnsCompartment(self):
        '''
        initialize specific columns for the SBtab type Compartment
        '''
        self.main_column    = None
        self.compart_column = None
        self.name_column    = None
        self.sbmlid_column  = None
        self.size_column    = None
        self.miriam_column  = None

        for i,column in enumerate(self.column_names):
            if column == '!Compartment':
                self.main_column    = i
                self.compart_column = i
            elif column == '!Name': self.name_column = i
            elif column == '!SBML::compartment::id': self.sbmlid_column = i
            elif column == '!Size': self.size_column = i
            elif column.startswith('!MiriamID'): self.miriam_column = i

    def initializeColumnsQuantityType(self):
        '''
        initialize specific columns for the SBtab type Compartment
        '''
        self.main_column = None
        self.qt_column   = None
        self.r_column    = None
        self.c_column    = None
        self.med_column  = None
        self.m_column    = None
        self.std_column  = None
        self.u_column    = None
        self.prov_column = None
        self.type_column = None
        self.src_column  = None
        self.logm_column = None
        self.logs_column = None
        self.min_column  = None
        self.max_column  = None
        self.pH_column   = None
        self.temp_column = None
        self.org_column  = None
        self.miriam_column  = None

        for i,column in enumerate(self.column_names):
            if column == 'QuantityType':
                self.main_column = i
                self.qt_column = i
            elif column == 'SBMLReactionID': self.r_column = i
            elif column == 'SBMLSpeciesID': self.c_column = i
            elif column == 'Median': self.med_column = i
            elif column == 'Mean': self.m_column = i
            elif column == 'Value': self.m_column = i
            elif column == 'Std': self.std_column = i
            elif column == 'Unit': self.u_column = i
            elif column == 'Provenance': self.prov_column = i
            elif column == 'Type': self.type_column = i
            elif column == 'Source': self.src_column = i
            elif column == 'logMean': self.logm_column = i
            elif column == 'logStd': self.logs_column = i
            elif column == 'Minimum': self.min_column = i
            elif column == 'Maximum': self.max_column = i
            elif column == 'pH': self.pH_column = i
            elif column == 'Temperature': self.temp_column = i
            elif column == 'OrganismName': self.org_name = i
            elif column.startswith('MiriamID'): self.miriam_column = i
            
    def checkColumns(self):
        '''
        checks whether the column names are known (according to the specification)
        '''
        specific_columns = self.table_type+'_columns.keys()'
        for column in self.column_names:
            if column in eval(specific_columns): pass
            else: print 'The column ',column,' is unknown according to the SBtab specification.'

