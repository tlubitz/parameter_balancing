#!/usr/bin/env python
import numpy
import scipy.linalg
import scipy.optimize
import copy
from . import bounded
import time

# declarations for the webversion defining some of the shown strings
name2index = {"standard chemical potential":0,
              "catalytic rate constant geometric mean":1,
              "Michaelis constant":2,
              "activation constant":3,
              "inhibitory constant":4,
              "concentration":5,
              "concentration of enzyme":6,
              "equilibrium constant":7,
              "substrate catalytic rate constant":8,
              "product catalytic rate constant":8,
              "forward maximal velocity":9,
              "reverse maximal velocity":10,
              "chemical potential":11,
              "reaction affinity":12}
index2displayname = {0:"Standard chemical potential (kJ/mol)",
                     1:"Catalytic rate constant geometric mean (1/s)",
                     2:"Michaelis constant (mM)",
                     3:"Activation constant (mM)",
                     4:"Inhibitory constant (mM)",
                     5:"Concentration (mM)",
                     6:"Concentration of enzyme (mM)",
                     7:"Equilibrium constant (dimensionless)",
                     8:"Catalytic rate constant (1/s)",
                     9:"Forward maximal velocity (mM/s)",
                     10:"Reverse maximal velocity (mM/s)",
                     11:"Chemical potential (kJ/mol)",
                     12:"Reaction affinity (kJ/mol)"}
header_names = ['QuantityType','SBMLReactionID','SBMLSpeciesID','Mode','Unit','Mean','Std','lnMean','lnStd']
inhibitory_sbos       = [20,206,207,536,537]
activation_sbos       = [13,21,459,461,462]
thermodynamic_indices = [0,11,12]
thermodynamics        = ['standard chemical potential','chemical potential','reaction affinity']

class ParameterBalancingError(Exception):
    '''
    default class for throwing individual pb errors
    '''
    def __init__(self, message):
        self.message = message
    def __str__(self):
        return self.message

class ParameterBalancing:
    '''
    class for the handling of parameter balancing
    '''
    def __init__(self,sbml_model,needed=True):
        '''
        initialise pb class
        '''
        self.model = sbml_model

        #initialise log file
        self.log   = 'Parameter balancing log file; %s\n\n'%(time.asctime())

        #is the uploaded model valid?
        if needed:
            try: self.model.getListOfSpecies()
            except: raise ParameterBalancingError('You have not used a valid SBML model.')
            self.gainModelInformation()

    def check4biomass(self):
        '''
        if the model appears to have a biomass function, we add a warning to the log file;
        biomass reactions should not have a reversible convenience kinetic and should be
        treated differently.
        '''
        warning  = False
        warning2 = False
        warning3 = False
        
        for species in self.model.getListOfSpecies():
            if 'biomass' in species.getName().lower() and not warning:
                self.log += '''### Warning about apparent biomass reaction: ###\nIt appears that the model has a biomass reaction (SBML ID %s). The convenience kinetics that are entered for the model reactions are not appropriate for biomass reactions. We suggest you replace the kinetics with something more appropriate like the one proposed by Hofmeyr et al. (2013). \n\n'''%(species.getId())
                warning = True                
        for reaction in self.model.getListOfReactions():
            if 'biomass' in reaction.getName().lower() and not warning2:
                self.log += '''### Warning about apparent biomass reaction: ###\nIt appears that the model has a biomass reaction (SBML ID %s). The convenience kinetics that are entered for the model reactions are not appropriate for biomass reactions. We suggest you replace the kinetics with something more appropriate like the one proposed by Hofmeyr et al. (2013). \n\n'''%(reaction.getId())
                warning2 = True
                continue
            if reaction.getNumReactants()>4 and not warning3:
                self.log += '''### Warning about apparent biomass reaction: ###\nThere is a reaction with an usually high amount of reactants (SBML ID %s). It may be a biomass reaction. The convenience kinetics that are entered for the model reactions are not appropriate for these reactions. We suggest you replace the kinetics with something more appropriate like the one proposed by Hofmeyr et al. (2013). \n\n'''%(reaction.getId())
                warning3 = True
                
    def gainModelInformation(self):
        '''
        this function gains initially all information about the SBML model:
        which species and reactions are contained? which modifiers?
        '''
        self.getParameterInformations()

        self.species_list    = []
        self.reaction_list   = []
        self.modifier_list   = []
        self.species2number  = {}
        self.reaction2number = {}
        self.medianed        = False      #decision variable, needed in several cases, so has to be declared very early
        
        #collect species and reactions from the SBML-model
        for i,species in enumerate(self.model.getListOfSpecies()):
            self.species_list.append(species.getId())
            self.species2number[species.getId()] = i+1

        for i,reaction in enumerate(self.model.getListOfReactions()):
            self.reaction_list.append(reaction.getId())
            for modifier in reaction.getListOfModifiers():
                self.modifier_list.append(modifier.getSpecies())
            self.reaction2number[reaction.getId()] = i+1
        
        (self.model_michaelis,self.model_inhibition,self.model_activation) = self.checkIfRequired()
        self.model_dict = {'Michaelis constant':self.model_michaelis,'inhibitory constant':self.model_inhibition,'activation constant':self.model_activation}
        self.buildReactionSpecifics()

        self.quantity2list = {'standard chemical potential':self.species_list,'chemical potential':self.species_list,'concentration':self.species_list,'catalytic rate constant geometric mean':self.reaction_list,'concentration of enzyme':self.reaction_list,'product catalytic rate constant':self.reaction_list,'substrate catalytic rate constant':self.reaction_list,'equilibrium constant':self.reaction_list,'forward maximal velocity':self.reaction_list,'reverse maximal velocity':self.reaction_list,'reaction affinity':self.reaction_list,'Michaelis constant':self.model_michaelis,'inhibitory constant':self.model_inhibition,'activation constant':self.model_activation}

        #add some information for the log file
        self.log += '### Model information ###\n'
        self.log += 'The SBML model has a total of %s Reactions and %s Species.\n'%(self.model.getNumReactions(),self.model.getNumSpecies())


    def getParameterInformations(self):
        '''
        this function reads a table file from the resources directory holding numerous informations on
        the handling of the parameters, the parameter details, and how to build the dependency matrix
        for the different parameter types
        '''
        #To make this module autonomous, I hardcode the usually dynamic file 'balancing_information'
        #open_file = open(libSBAnnotation.config.Config().getpath('resources','balancing_information'),'rw')
        #read_file = open_file.read()
        #header    = read_file.split('\n')[0].split('\t')
        #rows      = [row.split('\t') for row in read_file.split('\n')[1:]]

        header = ['QuantityType', 'Symbol', 'Unit', 'Constant', 'Related Element', 'Scaling', 'Dependence', 'Prior/pseudo mode', 'Prior/pseudo std.dev.', 'SBML element ', 'Abbreviation', 'MatrixInfo']
        rows = [['standard chemical potential', 'mu^0', 'kJ/mol', 'Thermodynamic', 'Species', 'Original', 'Basic', '-880', '4', 'Global parameter ', 'scp', '(1,0,0,0,0,0,0)'], ['catalytic rate constant geometric mean', 'k^V', '1/s', 'Kinetic', 'Reaction', 'Logarithmic', 'Basic', '10', '1', 'Local parameter', 'kcrg', '(0,1,0,0,0,0,0)'], ['Michaelis constant', 'k^M', 'mM', 'Kinetic', 'Reaction,Species', 'Logarithmic', 'Basic', '0.1', '1', 'Local parameter', 'kmc', '(0,0,1,0,0,0,0)'], ['activation constant', 'k^A', 'mM', 'Kinetic', 'Reaction,Species', 'Logarithmic', 'Basic', '0.1', '1', 'Local parameter', 'kac', '(0,0,0,1,0,0,0)'], ['inhibitory constant', 'k^I', 'mM', 'Kinetic', 'Reaction,Species', 'Logarithmic  ', 'Basic', '0.1', '1', 'Local parameter', 'kic', '(0,0,0,0,1,0,0)'], ['concentration of enzyme', 'u', 'mM', 'Dynamic', 'Reaction', 'Logarithmic', 'Basic', '0.00001', '1.5', 'Local parameter', 'eco', '(0,0,0,0,0,1,0)'], ['concentration', 'c', 'mM', 'Dynamic', 'Species', 'Logarithmic', 'Basic', '0.1', '1.5', 'Species (conc.)', 'c', '(0,0,0,0,0,0,1)'], ['equilibrium constant', 'k^eq', 'dimensionless', 'Thermodynamic', 'Reaction', 'Logarithmic', 'Derived', '1', '1.5', 'Local parameter', 'keq', '((-1/2.4942,"A"),0,0,0,0,0,0)'], ['substrate catalytic rate constant', 'k^cat+', '1/s', 'Kinetic', 'Reaction', 'Logarithmic', 'Derived', '10', '1.5', 'Local parameter  ', 'kcrf', '((-0.5/2.4942,"A"),1,(-0.5,"Z"),0,0,0,0)'], ['product catalytic rate constant', 'k^cat-', '1/s', 'Kinetic', 'Reaction', 'Logarithmic', 'Derived', '10', '1.5', 'Local parameter ', 'kcrr', '((0.5/2.4942,"A"),1,(0.5,"Z"),0,0,0,0)'], ['forward maximal velocity', 'v^max+', 'mM/s', 'Dynamic', 'Reaction', 'Logarithmic', 'Derived', '0.001', '2', 'Local parameter', 'vmaf', '((-0.5/2.4942,"A"),1,(-0.5,"Z"),0,0,1,0)'], ['reverse maximal velocity', 'v^max-', 'mM/s', 'Dynamic', 'Reaction', 'Logarithmic', 'Derived', '0.001', '2', 'Local parameter', 'vmar', '((0.5/2.4942,"A"),1,(0.5,"Z"),0,0,1,0)'], ['chemical potential', 'mu', 'kJ/mol', 'Dynamic', 'Species', 'Original', 'Derived', '-880', '680', '', '', "(1,0,0,0,0,0,(2.4942,'1'))"], ['reaction affinity', 'A', 'kJ/mol', 'Dynamic', 'Reaction', 'Original', 'Derived', '0', '10', '', '', '((-1,"A"),0,0,0,0,0,(-2.4942,"A"))'], ['']]
        
        for index,column_name in enumerate(header):
            if column_name == 'QuantityType': self.qt_column = index
            elif column_name == 'Unit': self.unit_column = index
            elif column_name == 'Constant': self.const_column = index
            elif column_name == 'Related Element': self.rel_column = index
            elif column_name == 'Scaling': self.scale_column = index
            elif column_name == 'Dependence': self.dep_column = index
            elif column_name == 'Prior/pseudo mode': self.median_column = index
            elif column_name == 'Prior/pseudo std.dev.': self.std_column = index
            elif column_name == 'Abbreviation': self.abbr_column = index
            elif column_name == 'MatrixInfo': self.matrix_column = index

        self.quantity2identifier      = {}
        self.quantity_type2unit       = {}
        self.quantity_type2mean_std   = {}
        self.quantity_type2median_std = {}
        self.prior_list          = []
        self.pseudo_list         = []
        self.quantity_types      = []
        self.species_parameters  = []
        self.reaction_parameters = []
        self.reaction_species_parameters  = []
        self.condition_dependent_quantity = []
        
        for row in rows:
            if len(row) == len(header):
                self.quantity_types.append(row[0])
                self.quantity2identifier[row[0]] = row[self.abbr_column]
                self.quantity_type2unit[row[0]] = row[self.unit_column]
                if row[self.unit_column] == 'kJ/mol':
                    self.quantity_type2mean_std[row[0]] = [float(row[self.median_column]),float(row[self.std_column])]
                else: self.quantity_type2median_std[row[0]] = [float(row[self.median_column]),float(row[self.std_column])]
                if row[self.dep_column] == 'Basic': self.prior_list.append(row[self.qt_column])
                elif row[self.dep_column] == 'Derived': self.pseudo_list.append(row[self.qt_column])
                if row[self.rel_column] == 'Species': self.species_parameters.append(row[self.qt_column])
                elif row[self.rel_column] == 'Reaction':
                    self.reaction_parameters.append(row[self.qt_column])
                elif row[self.rel_column] == 'Reaction,Species':
                    self.reaction_species_parameters.append(row[self.qt_column])
                if row[self.const_column] == 'Thermodynamic': self.condition_dependent_quantity.append(row[self.qt_column])
            
    def checkIfRequired(self):
        '''
        not every reaction/species-combo needs an activation, inhibitory, or Michaelis constant.
        this function checks, which combos DO need them
        '''
        Michaelis  = []
        inhibitors = []
        activators = []
        
        for reaction in self.model.getListOfReactions():
            #first: check for inhibitory and activation constants (only identifiable by SBO terms)
            #Update: no inhibitory constants suck. So forget about SBO terms and assign one to each reaction, as initially thought
            #Thus -> required or not, I add both to "Reaction"-classifier in information.tsv
            for modifier in reaction.getListOfModifiers():
                if modifier.getSBOTerm() in inhibitory_sbos: # and modifier.getSpecies() in self.modifier_list:
                    inhibitors.append(['inhibitory constant',reaction.getId(),modifier.getSpecies()])
                elif modifier.getSBOTerm() in activation_sbos: # and species.getId() in self.modifier_list:
                    activators.append(['activation constant',reaction.getId(),modifier.getSpecies()])
            '''
            for species in self.model.getListOfSpecies():
                if species.getSBOTerm() in inhibitory_sbos and species.getId() in self.modifier_list:
                    inhibitors.append(['inhibitory constant',reaction.getId(),species.getId()])
                elif species.getSBOTerm() in activation_sbos and species.getId() in self.modifier_list:
                    activators.append(['activation constant',reaction.getId(),species.getId()])

            '''
            #second: check for Michaelis constants
            for reactant in reaction.getListOfReactants():
                Michaelis.append(['Michaelis constant',reaction.getId(),reactant.getSpecies()])
            for product in reaction.getListOfProducts():
                Michaelis.append(['Michaelis constant',reaction.getId(),product.getSpecies()])

        self.model_specific = Michaelis+inhibitors+activators
        return Michaelis,inhibitors,activators

    def buildReactionSpecifics(self):
        '''
        generates some dictionaries linking the reactions with its reactants, products and stoichiometric coefficients
        '''
        self.reactions_reactants = {}
        self.reactions_products  = {}
        
        #generate some dictionaries to link the reactions to the reactants/products
        #and their stoichiometric coefficients
        for i,reaction in enumerate(self.model.getListOfReactions()):
            reactants = []
            stoich    = []
            for reactant in reaction.getListOfReactants():
                reactants.append(reactant.getSpecies())
                thisstoich = reactant.getStoichiometry()
                if thisstoich != thisstoich:
                    thisstoich = 1.
                stoich.append(thisstoich*(-1))
            self.reactions_reactants[reaction.getId()] = (reactants,stoich)

            products = []
            stoich   = []
            for product in reaction.getListOfProducts():
                products.append(product.getSpecies())
                thisstoich = product.getStoichiometry()
                if thisstoich != thisstoich:
                    thisstoich = 1.0
                stoich.append(thisstoich)
            self.reactions_products[reaction.getId()] = (products,stoich)

    def makeEmptySBtab(self,pmin,pmax,parameter_dict):
        '''
        if no SBtab is given, create an empty one for the model using
        function makeSBtab and some default parameters
        '''
        sbtab_string = 'QuantityType\tSBMLReactionID\tSBMLSpeciesID\tMode\tUnit\tMean\tStd\n\t\t\t\t\t\t\n'
        header = header_names
        filename = 'newSBtab.tsv'
        organism = 'All organisms'
        volume   = 42.0
        (new_header,rows) = self.makeSBtab(sbtab_string,header,filename,organism,volume,pmin,pmax,parameter_dict)

        return new_header,rows


    def makeSBtab(self,SBtabFile,header,fileName,organism,volume,pmin,pmax,parameter_dict):
        '''
        makes an insertable SBtab-file out of a possibly erraneous SBtab-File
        '''
        #if the new SBtabFormat is given, update the variables
        self.header       = header
        inconsistent_rows = self.makeStringToLists(SBtabFile,fileName)

        self.organism     = organism
        self.volume       = volume
        self.new_header   = header_names
        self.new_rows     = []

        self.pmin         = pmin
        self.pmax         = pmax
        self.parameter_dict = parameter_dict
        
        self.getColumnIndices(self.header)
        self.rows = self.tidyUpSBtab(inconsistent_rows,self.organism,False)

        for header in header_names:
            if not header in self.header:
                self.header.append(header)
                for row in self.rows:
                    row.append('')

        checklength = len(self.header)       
        self.available_parameters = []
        for i,row in enumerate(self.rows):
            if len(row) == checklength:
                if row[self.qt_column] in self.species_parameters:
                    self.available_parameters.append([row[self.qt_column],row[self.s_column]])
                elif row[self.qt_column] in self.reaction_parameters:
                    self.available_parameters.append([row[self.qt_column],row[self.r_column]])
                elif row[self.qt_column] in self.reaction_species_parameters:
                    self.available_parameters.append([row[self.qt_column],row[self.r_column],row[self.s_column]])

        #START GENERATING THE NEW SBtab-FILE:
        for species in self.species_list:
            for quantity in self.species_parameters:
                if [quantity,species] in self.available_parameters:
                    self.new_rows.append(self.existingRow(species,quantity))
                else: self.new_rows.append(self.newRow(species,quantity))
                
        for reaction in self.reaction_list:
            for quantity in self.reaction_parameters:
                if [quantity,reaction] in self.available_parameters:
                    self.new_rows.append(self.existingRow(reaction,quantity))
                else:
                    self.new_rows.append(self.newRow(reaction,quantity))        
                    
        tuple_list = [self.model_michaelis,self.model_activation,self.model_inhibition]
        for i,liste in enumerate(tuple_list):
            for reaction_species in liste:
                if reaction_species in self.available_parameters:
                    self.new_rows.append(self.existingRow(reaction_species,self.reaction_species_parameters[i]))
                else:
                    self.new_rows.append(self.newRow(reaction_species,self.reaction_species_parameters[i]))

        uniques = []
        for entry in self.available_parameters:
            if not entry in uniques: uniques.append(entry)

        self.log += 'You have used %s values to describe %s unique kinetic parameters as input for the balancing.\n'%(len(self.rows),len(self.new_rows))
        self.log += 'The parameter balancing will output %s standard chemical potentials, %s catalytic rate constants geometric mean, %s activation constants, %s inhibitory constants, %s concentrations, %s concentrations of enzyme, %s equilibrium constants, %s substrate catalytic rate constants, %s product catalytic rate constants, %s forward maximal velocities, %s reverse maximal velocities, %s chemical potentials, and %s reaction affinities.\n'%(self.model.getNumSpecies(),self.model.getNumReactions(),self.model.getNumReactions(),self.model.getNumReactions(),self.model.getNumSpecies(),self.model.getNumReactions(),self.model.getNumReactions(),self.model.getNumReactions(),self.model.getNumReactions(),self.model.getNumReactions(),self.model.getNumReactions(),self.model.getNumSpecies(),self.model.getNumReactions())

        return (self.new_header,self.new_rows)

    def getColumnIndices(self,header):
        '''
        this function collects the indices of the different columns
        '''
        self.temp_column = False
        self.pH_column   = False
        self.src_column  = False
        self.min_column  = False
        self.max_column  = False

        for i,name in enumerate(header):
            if name == 'Quantity': self.q_column = i
            elif name == 'QuantityType': self.qt_column = i
            elif name == 'SBMLReactionID': self.r_column = i               
            elif name == 'SBMLSpeciesID': self.s_column = i
            elif name == 'Mean': self.m_column = i
            elif name == 'Value': self.m_column = i
            elif name == 'Unit': self.u_column = i
            elif name == 'Std': self.std_column = i
            elif name == 'Mode': self.med_column = i                
            elif name == 'Temperature': self.temp_column = i
            elif name == 'pH': self.pH_column = i
            elif name == 'Source' or name == 'Reference': self.src_column = i
            elif name == 'OrganismName': self.o_column = i
            elif name == 'Type': self.t_column = i
            elif name == 'lnMean': self.logm_column = i
            elif name == 'lnStd': self.logs_column = i
            elif name == 'Minimum': self.min_column = i
            elif name == 'Min': self.min_column = i
            elif name == 'Maximum': self.max_column = i
            elif name == 'Max': self.max_column = i

        if self.temp_column != False and not 'Temperature' in self.new_header:
            self.new_header.append('Temperature')
        if self.pH_column != False and not 'pH' in self.new_header:
            self.new_header.append('pH')
        if self.min_column != False and not 'Min' in self.new_header:
            self.new_header.append('Min')
        if self.max_column != False and not 'Max' in self.new_header:
            self.new_header.append('Max')

    def makeStringToLists(self,SBtabFile,fileFormat):
        '''
        this function makes the badly editable tsv-string to a lists-tuple.
        every tuple is one row.
        '''
        all_rows = SBtabFile.split('\n')
        rows     = []
        if fileFormat.endswith('tsv'):
            for row in all_rows[1:]:
                rows.append(row.split('\t'))
        elif fileFormat.endswith('csv'):
            for row in all_rows[1:]:
                rows.append(row.split('\t'))
        else:
            raise ParameterBalancingError('The given fileformat for the SBtabfile is not supported. Please insert a .tsv or .csv-file instead of ',fileFormat)

        return rows
    
    def tidyUpSBtab(self,inconsistent_rows,organism,ig_bounds):
        '''
        removes a lot of stuff from SBtab that is either cramping the functionality or
        deletes entries that must not be used in parameter balancing
        '''
        consistent_rows   = []                        #rows that fulfill all consistency criteria
        log_header        = False

        for row in inconsistent_rows:
            if len(row) == len(inconsistent_rows[0]) and ''.join(row) != '':
                #1st: if the user has specified an organism then remove all others
                if organism and organism != 'All organisms':
                    if row[self.o_column] != self.organism or row[self.o_column] != '':
                        continue
                                            
                #2nd: we do not use KMedDB entries anymore!
                if row[self.src_column] == 'KMedDB':
                    continue

                #3rd: inhibitory constants need reaction AND species! Same for Michaelis constants. (Update: only MCs)
                #if row[self.qt_column] == 'inhibitory constant' or row[self.qt_column] == 'Michaelis constant' or row[self.qt_column] == 'activation constant':
                if row[self.qt_column] == 'Michaelis constant':
                    if row[self.r_column] not in self.reaction_list or row[self.s_column] not in self.species_list:
                        continue

                #4th: if we are dealing with the unit "molecules/cell" on enzyme concentrations, recalculate the unit
                if row[self.qt_column] == 'concentration of enzyme' and row[self.u_column] == 'molecules/cell':
                    value = float(row[self.m_column])
                    row[self.m_column] = str(value*0.00000004)
                    row[self.std_column] = str(value*0.00000002)
                    row[self.u_column] = 'mM'

                #5th: we do not always want to ask for the content of empty columns. so we
                #agree to set all of them to an empty string
                for i,value in enumerate(row):
                    if value == 'nan' or value == '-' or value == 'NaN' or value == 'None' or value == None:
                        row[i] = ''

                #5bth: Exclude values that lie outside of the given boundaries, if requested
                if 'boundary_values' in self.parameter_dict and ig_bounds:
                    if self.parameter_dict['boundary_values'] == 'ignore':
                        if row[self.qt_column] != '' and self.pmin[row[self.qt_column]] != None:
                            if float(row[self.m_column]) < float(self.pmin[row[self.qt_column]]):
                                if log_header == False:
                                    self.log += '\n### Warnings about ignored values that lie out of the boundaries ###\n'
                                    log_header = True
                                self.log += 'The value %s for the %s of %s,%s lies below the requested minimum value of %s. It is ignored for the balancing.\n'%(row[self.m_column],row[self.qt_column],row[self.r_column],row[self.s_column],self.pmin[row[self.qt_column]])
                                continue
                        if row[self.qt_column] != '' and self.pmax[row[self.qt_column]] != None:
                            if float(row[self.m_column]) > float(self.pmax[row[self.qt_column]]):
                                if log_header == False:
                                    self.log += '\n### Warnings about ignored values that lie out of the boundaries ###\n'
                                    log_header = True
                                self.log += 'The value %s for the %s of %s,%s lies above the requested maximum value of %s. It is ignored for the balancing.\n'%(row[self.m_column],row[self.qt_column],row[self.r_column],row[self.s_column],self.pmax[row[self.qt_column]])
                                continue
                        
                #some specific stuff for the balancing with fixed eq constants and concentrations
                '''
                if row[self.qt_column] == 'concentration':
                    row[self.min_column] = 0.001
                    row[self.max_column] = 100
                    row[self.std_column] = 1
                if row[self.qt_column] == 'equilibrium constant':
                    row[self.std_column] = 1
                if row[self.qt_column] == 'reaction affinity':
                    row[self.std_column] = 0.01
                if row[self.qt_column] == 'standard chemical potential':
                    row[self.std_column] = 1
                
                if row[self.qt_column] == 'equilibrium constant' or row[self.qt_column] == 'concentration':
                    row[self.std_column] = float(0.00000001)
                    row[self.min_column] = row[self.m_column]
                    row[self.max_column] = row[self.m_column]
                else:
                    row[self.min_column] = ''
                    row[self.max_column] = ''
                '''
                if row[self.m_column] == '':
                    continue
                                  
                #6th: it leads to followup-mistakes if reaction-parameters also have species specified or vice versa
                #-> so we erase this eventuality
                if row[self.qt_column] == 'standard chemical potential' or row[self.qt_column] == 'concentration' or row[self.qt_column] == 'chemical potential':
                    row[self.r_column] = ''
                elif row[self.qt_column] == 'catalytic rate constant geometric mean' or row[self.qt_column] == 'concentration of enzyme' or row[self.qt_column] == 'equilibrium constant' or row[self.qt_column] == 'substrate catalytic rate constant' or row[self.qt_column] == 'product catalytic rate constant' or row[self.qt_column] == 'forward maximal velocity' or row[self.qt_column] == 'reverse maximal velocity' or row[self.qt_column] == 'reaction affinity':
                    row[self.s_column] = ''

                if row[self.std_column] == '':
                    row[self.std_column] = float(row[self.m_column])*0.5

                #last: if everything's alright, append the row to the set of consistent rows
                consistent_rows.append(row)
                
        return consistent_rows

    def newRow(self,name,quantity):
        '''
        generates one reaction row that is required
        '''
        row = ['']*len(self.new_header)
        row[0] = quantity
        if quantity in self.reaction_species_parameters:
            row[1] = name[1]
            row[2] = name[2]
        elif quantity in self.species_parameters: row[2] = name
        else: row[1] = name
        row[4] = self.quantity_type2unit[quantity]

        return row

    def existingRow(self,name,quantity):
        '''
        generates one species row that is required and available
        '''
        #first: check whether there are multiple entries for this specific parameter. if so, do mean
        #(if we have a quantity with measurement conditions, we must not mean them)
        value_dict = False
        rands      = False
        amount     = 1
        if quantity in self.reaction_species_parameters: rands = True
        if rands:                #do we need reaction AND species?
            if self.available_parameters.count([quantity,name[1],name[2]])>1:
                value_dict = self.meanRow(name,quantity)
        elif self.available_parameters.count([quantity,name])>1:   #elsewise
            value_dict = self.meanRow(name,quantity)
        
        #next: build the parameter row
        used_rows = []
        new_row = ['']*len(self.new_header)

        for row in self.rows:
            for i in range(0,amount):
                further = False
                if rands:
                    if row[self.qt_column] == quantity and row[self.s_column] == name[2] and row[self.r_column] == name[1]:
                        new_row[1] = row[self.r_column]
                        new_row[2] = row[self.s_column]
                        further    = True
                else:
                    if row[self.qt_column] == quantity and (row[self.s_column] == name or row[self.r_column] == name) and not row in used_rows:
                        new_row[1] = row[self.r_column]
                        new_row[2] = row[self.s_column]
                        further    = True
                if further:
                    new_row[0] = row[self.qt_column]
                    if value_dict:
                        new_row[3] = str(value_dict['Mode'])
                        new_row[5] = str(value_dict['Mean'])
                        new_row[6] = str(value_dict['Std'])
                    else:
                        new_row[5] = row[self.m_column]
                        new_row[6] = row[self.std_column]
                        new_row[3] = str(round(numpy.exp(self.normal2log([float(new_row[5])],[float(new_row[6])],[new_row[0]])[0])[0],4))
                    if quantity in self.quantity_type2mean_std.keys(): new_row[3] = new_row[5]
                    new_row[4] = self.quantity_type2unit[row[self.qt_column]]
                    
                    #columns seven and eight are the logmean and logstd
                    #optional columns (columns 9 and more)
                    j = 9
                    if self.temp_column:
                        new_row[j] = row[self.temp_column]
                        j += 1
                    if self.pH_column:
                        new_row[j] = row[self.pH_column]
                        j += 1
                    if self.min_column and self.max_column:
                        new_row[j] = row[self.min_column]
                        new_row[j+1] = row[self.max_column]
                         
                    if amount>1 and i<amount:
                        multi_rows.append(new_row)
                        new_row = ['']*len(self.new_header)
                        used_rows.append(row)
                        if len(multi_rows)==amount:
                            return multi_rows

        return new_row

    def meanRow(self,name,quantity):
        '''
        if there are more than one value for one parameter, calculate the mean
        '''
        #collect available means and stds
        means      = []
        stds       = []
        rs         = False
        if quantity in self.reaction_species_parameters: rs = True
        for row in self.rows:
            if rs:
                if row[self.qt_column] == quantity and row[self.r_column] == name[1] and row[self.s_column] == name[2]:
                    means.append(float(row[self.m_column]))
                    stds.append(float(row[self.std_column]))
            else:
                if row[self.qt_column] == quantity and (row[self.s_column] == name or row[self.r_column] == name):
                    means.append(float(row[self.m_column]))
                    stds.append(float(row[self.std_column]))
        #build the mean
        #######
        #we have a problem here: actually, for thermodynamic quantities, we want to use the geometric mean instead of the
        #arithmetic mean. But, as it turns out, the .gmean() function of scipy employs the .log() function, causing a crash
        #when faced with values lower than zero (which are common in thermodynamic quantities such as Gibbs energies.
        #Thusly, I am trying to use the arithmetic mean for all quantities now, also the thermodynamic ones.
        #######
        
        #if quantity in self.quantity_type2mean_std.keys():
        #    #build geometric mean
        #    import scipy.stats
        #    denominator = 0
        #    for std in stds:
        #        if std != 0.0: denominator = denominator + (1/std**2)
        #        else: denominator = denominator + (1/numpy.log(2)**2)
        #
        #    if denominator == 0: denominator += 0.000000001
        #    mean   = scipy.stats.gmean(means)
        #    std    = numpy.exp(numpy.sqrt(1/denominator))
        #    if not quantity in self.quantity_type2mean_std.keys(): median = numpy.exp(self.normal2log([mean],[std],False)[0])[0]
        #    else: median = mean
        #    value_dict = dict([('Mean',mean),('Std',std),('Median',median)])
        
        #build arithmetic mean
        denominator = 0
        for i,std in enumerate(stds):
            if std != 0.0: denominator = denominator + (1/std**2)
            else: denominator = denominator + (1/numpy.log(2)**2)
        if denominator == 0: denominator += 0.000000001
        std = numpy.sqrt(1/denominator)

        numerator   = 0
        denominator = 0
        for i,mean in enumerate(means):
            numerator   = numerator + mean/stds[i]**2
            denominator = denominator + 1/stds[i]**2
        if denominator == 0: denominator += 0.000000001
        mean   = numerator/denominator
        if not quantity in self.quantity_type2mean_std.keys(): median = numpy.exp(self.normal2log([mean],[std],False)[0])[0]
        else: median = mean            
        value_dict = dict([('Mean',mean),('Std',std),('Mode',median)])

        return value_dict

    def fillSBtab(self,sbtab_file,priors=None,pseudos=None):
        '''
        fills the values in the given SBtabfile
        '''
        self.prior_used  = False
        self.pseudo_used = False

        first_row = sbtab_file[0]
        if priors:
            self.prior_used = True
            for row in sbtab_file:
                if len(row) == len(first_row):
                    if row[5] == '':
                        try:
                            row[5] = str(priors[row[0]][0])
                            row[6] = str(priors[row[0]][1])
                        except: pass

        if pseudos:
            self.pseudo_used = True
            for row in sbtab_file:
                if len(row) == len(first_row):
                    if row[5] == '':
                        try:
                            row[5] = str(pseudos[row[0]][0])
                            row[6] = str(pseudos[row[0]][1])
                        except: pass



        (self.log_means,self.log_stds) = self.makeDefaultTable()

        #generate a new prior
        if priors:
            means = [None]*14
            stds  = [None]*14
            means[0] = priors['standard chemical potential'][0]
            stds[0] = priors['standard chemical potential'][1]
            means[1] = priors['catalytic rate constant geometric mean'][0]
            stds[1] = priors[ 'catalytic rate constant geometric mean'][1]
            means[2] = priors['Michaelis constant'][0]
            stds[2] = priors['Michaelis constant'][1]
            means[3] = priors['activation constant'][0]
            stds[3] = priors['activation constant'][1]       
            means[4] = priors['inhibitory constant'][0]
            stds[4] = priors['inhibitory constant'][1]
            means[5] = priors['concentration'][0]
            stds[5] = priors['concentration'][1]
            means[6] = priors['concentration of enzyme'][0]
            stds[6] = priors['concentration of enzyme'][1]
            if pseudos:
                for [name, index] in [['equilibrium constant', 7],
                                      ['substrate catalytic rate constant', 8],
                                      ['product catalytic rate constant', 9],
                                      ['forward maximal velocity', 10],
                                      ['reverse maximal velocity', 11],
                                      ['chemical potential', 12],
                                      ['reaction affinity', 13]]:
                    if name in pseudos:
                        means[index] = pseudos[name][0]
                        stds[index] = pseudos[name][1]
            else:
                for [name, index] in [['equilibrium constant', 7],
                                      ['substrate catalytic rate constant', 8],
                                      ['product catalytic rate constant', 9],
                                      ['forward maximal velocity', 10],
                                      ['reverse maximal velocity', 11],
                                      ['chemical potential', 12],
                                      ['reaction affinity', 13]]:
                    if name in self.quantity_type2median_std.keys():
                        (mean,std) = self.med10std2normal([self.quantity_type2median_std[name][0]],[self.quantity_type2median_std[name][1]],[name])
                        means[index] = mean[0]
                        stds[index]  = stds[0]
                    else:
                        (mean,std) = self.med10std2normal([self.quantity_type2mean_std[name][0]],[self.quantity_type2mean_std[name][1]],[name])
                        means[index] = mean[0]
                        stds[index]  = stds[0]
            types = self.prior_list+self.pseudo_list

            for i in range(means.count(None)):
                means.remove(None)

            for i in range(stds.count(None)):
                stds.remove(None)
                
            (self.means,self.log_stds) = self.normal2log(means,stds,types)

        return sbtab_file


    def calcMedians(self,header,rows):
        '''
        if only the default values are entered into the SBtab, we still have to calculate the medians
        for every row. This is done by this function.
        '''
        for i,column in enumerate(header):
            if column == 'QuantityType': qt_column = i
            if column == 'Mean'  : mean_column = i
            if column == 'Std'   : std_column  = i
            if column == 'Mode': med_column  = i
        
        for row in rows:
            if row[qt_column] in self.quantity_type2mean_std.keys():
                row[med_column] = row[mean_column]
            else:
                (log_mean,log_std) = self.normal2log([float(row[mean_column])],[float(row[std_column])],row[qt_column])
                row[med_column] = numpy.exp(log_mean[0])
        return rows
        

    def makeDefaultTable(self):
        '''
        generate a table of values for every parameter for every reaction/species
        '''
        medians    = []
        std10logs  = []
        quantities = []
        types      = []

        for parameter_type in self.quantity_type2mean_std.keys():
            medians.append(self.quantity_type2mean_std[parameter_type][0])
            quantities.append(parameter_type)
            std10logs.append(self.quantity_type2mean_std[parameter_type][1])
            types.append(parameter_type)
        for parameter_type in self.quantity_type2median_std.keys():
            medians.append(self.quantity_type2median_std[parameter_type][0])
            std10logs.append(self.quantity_type2median_std[parameter_type][1])
            quantities.append(parameter_type)
            types.append(parameter_type)       
            
        #make the median/std10log values to normal  values
        (means,stds) = self.med10std2normal(medians,std10logs,types)

        self.prior_values = {}
        for i,quantity in enumerate(quantities):
            self.prior_values[quantity] = [(means[i],stds[i])]

        return means,stds

    def normal2log(self,means,stds,types):
        '''
        generates log values for normal values
        '''
        log_means = []
        log_stds  = []

        for i,mean in enumerate(means):
            if types and types[i] in self.quantity_type2mean_std.keys():
                if mean == 0.0: mean += 0.00001
                if stds[i] == 0.0: stds[i] += 0.0001
                log_means.append(float(mean))
                log_stds.append(float(stds[i]))
            else:
                if mean == 0.0: mean += 0.00001
                if stds[i] == 0.0: stds[i] += 0.00001
                term = numpy.log(1+numpy.square(float(stds[i])/float(mean)))
                if term == 0.0: term += 0.00001
                numpy.seterr(invalid='raise')
                log_means.append(numpy.log(float(mean))-0.5*term)
                log_stds.append(numpy.sqrt(numpy.log(1+numpy.square(float(stds[i])/float(mean)))))
        if 'nan' in log_means:
            raise ParameterBalancingError('The logarithm of one of your given mean values is invalid.')

        return log_means,log_stds

    def printWarning(self,type,flag):
        print('There was an %s error in the numerics. This may be due to broadly chosen probability distributions. Try the usage of pseudo values for a fix.'%(type))  

    def log2normal(self,log_means,log_stds,types=None):
        '''
        generates a list of the normal values from a list of the log values
        '''
        means = []
        stds  = []

        numpy.seterrcall(self.printWarning)
        numpy.seterr(all='call')
        for i,log_mean in enumerate(log_means):
            if types and types[i] in self.quantity_type2mean_std.keys():
                if log_mean == 0.0: log_mean += 0.000000001
                if log_stds[i] == 0.0: log_stds[i] += 0.00000001
                means.append(log_mean)
                stds.append(log_stds[i])
            else:
                if log_mean == 0.0: log_mean += 0.00000001
                if log_stds[i] == 0.0: log_stds[i] += 0.00000001
                means.append(numpy.exp(log_mean+0.5*numpy.square(log_stds[i])))
                stds.append(numpy.sqrt((numpy.exp(numpy.square(float(log_stds[i])))-1)*(numpy.exp(2*log_mean+numpy.square(float(log_stds[i]))))))
                #try: stds.append(numpy.sqrt((numpy.exp(numpy.square(float(log_stds[i])))-1)*(numpy.exp(2*log_mean+numpy.square(float(log_stds[i]))))))
                #except: stds.append(None)
        return means,stds

    def med10std2normal(self,medians,stdlogs,types):
        '''
        this is the users default choice for missing values: inserted are the median and
        the stdlog10, the output are the corresponding mean value and the standard dev.
        '''
        #first, we make the median to the mean_ln and the corresponding std_ln
        log_means = []
        log_stds  = []

        for i,median in enumerate(medians):
            if types and types[i] in self.quantity_type2mean_std.keys():
                log_means.append(median)
                log_stds.append(stdlogs[i])
            else:
                log_means.append(numpy.log(median))
                log_stds.append(stdlogs[i]/(1/numpy.log(10)))

        #then, we use the two log-values to reconstruct the mean and std:
        (means,stds) = self.log2normal(log_means,log_stds,types)

        return means,stds

    def addConfig2log(self):
        '''
        adds config information for the balancing
        '''
        self.log += '\n### You have used an options file for the parameter balancing. ###\nThe options are as follows:\n'
        self.log += '!Option\t!Value\n'
        for entry in self.parameter_dict.keys():
            self.log += '%s\t%s\n'%(entry,self.parameter_dict[entry])

    def makeBalancing(self,parameter_dict,new_sbtab,old_sbtab,pmin,pmax):
        '''
        generates the values for the parameter balancing
        '''
        #initialise log file
        self.parameter_dict = parameter_dict
        if 'config' in self.parameter_dict:
            self.addConfig2log()
        
        #initialize needed variables
        self.old_header       = old_sbtab.split('\n')[0].split('\t')
        self.pmin             = pmin
        self.pmax             = pmax
        inconsistent_old_rows = self.makeStringToLists(old_sbtab,'bla.tsv')
        self.new_header       = header_names
        self.new_rows         = new_sbtab
        self.sheet            = self.getSheet()
        self.temperature      = float(self.parameter_dict['Temperature'])
        self.pH               = float(self.parameter_dict['pH'])

        #build needed vectors and matrices
        self.desired_parameters = self.buildDesiredParameters()
        self.x_vector           = self.collectAvailableValues(inconsistent_old_rows)
        self.theta_vector       = self.buildThetaVector()
        (self.C_0,self.C_x)     = self.buildCovarianceMatrices()


        #build the dependence matrix D and the data specific D_x
        self.D   = self.buildDependenceMatrix()
        self.D_x = self.buildSpecificDependenceMatrix()

        #KEY PART: calculating the posteriori covariance matrix C_post and the posteriori mean vector mean_post
        self.calculatePosteriori()

        #from here: post-stuff, won't be a problem
        #extract the stds from the covariance matrix C_post for export into GUI
        self.log_stds_post = self.extractPostCovariance()

        #make normal values again
        (self.mean_post,self.stds_post) = self.log2normal(self.log_mean_post,self.log_stds_post,self.type_order)
        ################################################################
        #generating minimization problem
        self.old_means = copy.deepcopy(self.mean_post)
        
        if self.bounds_inc.count((None,None)) != len(self.bounds_inc):# and False:                #and False ---> mustn't be loaded in web interface
            #generating the medians to bound them                                     #because the files have to be generated manually!!
            #print('in optimiser')                                                      #only for super-cool users...
            medians = []
            medstds = []

            (self.means_inc,self.stds_inc) = self.log2normal(self.means_log_inc,self.stds_log_inc,self.types_inc)

            for i,value in enumerate(self.means_inc):
                if not self.types_inc[i] in self.quantity_type2mean_std.keys():
                    (log_mean,log_std) = self.normal2log([self.means_inc[i]],[self.stds_inc[i]],self.types_inc[i])
                    medians.append(numpy.exp(log_mean[0]))
                else:
                    medians.append(self.means_inc[i])
                medstds.append(max(self.stds_inc[i],10))

            #optimization function
            def log_mean_post_func(q):
                return numpy.dot(((q-medians).transpose()),(numpy.dot(numpy.linalg.inv(self.C_post_inc),(q-medians))))

            #setting proper boundaries for parameters that have no boundaries set in SBtab
            new_boundaries = []
            is_logarithmic = []

            for i,bound in enumerate(self.bounds_inc):
                if bound == ('','') or bound == (None,None):
                    #setting boundaries for thermodynamic parameters (kJ/mol)
                    if name2index[self.types_inc[i]] in thermodynamic_indices:
                        is_logarithmic.append(False)
                        if (max(medians[i]-medstds[i]*2,-3000))<(min(3000,medians[i]+medstds[i]*2)):
                            new_boundaries.append((max(medians[i]-medstds[i]*2,-3000),min(3000,medians[i]+medstds[i]*2)))
                        else:
                            new_boundaries.append((medians[i]+medstds[i]*2,medians[i]-medstds[i]*2))
                    #setting boundaries for all other parameters
                    else:
                        is_logarithmic.append(True)
                        if (medians[i]-medstds[i]*4)<(medians[i]+medstds[i]*4):
                            new_boundaries.append((max(medians[i]-medstds[i]*4,0.00001),medians[i]+medstds[i]*4))
                        else:
                            new_boundaries.append((medians[i]+medstds[i]*4,medians[i]-medstds[i]*4))
                else:
                    is_logarithmic.append(True)
                    new_boundaries.append((float(bound[0]),float(bound[1])))

            proper_boundaries = []

            for boundaries in new_boundaries:
                new_bound = []
                if boundaries[0] == float(0.0):
                    new_bound.append(0.00001)
                else:
                    new_bound.append(boundaries[0])
                if boundaries[1] == float(0.0):
                    new_bound.append(0.00001)
                else:
                    new_bound.append(boundaries[1])
                proper_boundaries.append(new_bound)

           
            f = open('medians.txt','w')
            for i,element in enumerate(medians):
                if not i == len(medians)-1:
                    f.write(str(element)+',')
                else:
                    f.write(str(element))
            f.close()

            g = open('cpost.txt','w')
            for line in self.C_post_inc:
                for i,element in enumerate(line):
                    if not i == len(line)-1:
                        g.write(str(element)+',')
                    else:
                        g.write(str(element))
                g.write('\n')
            g.close()

            #generating optimization and updating responsible mean vector
            new_medians = bounded.fmin_gen(log_mean_post_func,numpy.array(medians),population_size=10,survivors=3,generations=5,bounds=proper_boundaries,use_pp=False,variable_is_logarithmic=is_logarithmic,disp=0)
            #new_medians = libTI2.optimize.bounded.fmin_differential_evolution(log_mean_post_func,numpy.array(medians),population_size=100,generations=100,bounds=proper_boundaries,variable_is_logarithmic=is_logarithmic,disp=1)

            for i,value in enumerate(new_medians):
                if value+0.001 < proper_boundaries[i][0] or value-0.001 > proper_boundaries[i][1]:
                    self.types_inc[i]
                    print('NEW_MODE value out of bound: ',value,' [bounds: ',proper_boundaries[i],']')

            (new_medians_log,new_stds_log) = self.normal2log(new_medians,self.stds_inc,self.types_inc)

            self.C_post = numpy.dot((numpy.dot(self.D,self.C_post_inc)),self.D.transpose())
            self.log_mean_post = numpy.dot(self.D,new_medians_log)
            self.log_stds_post = self.extractPostCovariance()

            (self.mean_post_opt,self.stds_post_opt) = self.log2normal(self.log_mean_post,self.log_stds_post,self.type_order)
            
            self.medianed  = True
            #self.mean_post = copy.deepcopy(self.mean_post_a)
      
        #################################################################
        #make value-dictionaries to realize the insertion of the computed values into the SBtab-file in GUI
        self.balancedSBtab = self.getNewSBtab()
        C_string = self.makeCpostString()
        shannons = self.getShannons()

        return self.mean_post,self.mean_post_inc,C_string,self.C_post_inc,self.D,shannons,self.balancedSBtab,self.log

    def getSheet(self):
        '''
        get the sheet that tells us, how to build up which matrix and further stuff
        '''
        sheet = {"equilibrium constant":((-1/2.4942,"A"),0,0,0,0,0,0), "substrate catalytic rate constant":((-0.5/2.4942,"A"),1,(-0.5,"Z"),0,0,0,0), "product catalytic rate constant":((0.5/2.4942,"A"),1,(0.5,"Z"),0,0,0,0), "forward maximal velocity":((-0.5/2.4942,"A"),1,(-0.5,"Z"),0,0,1,0), "reverse maximal velocity":((0.5/2.4942,"A"),1,(0.5,"Z"),0,0,1,0),"chemical potential":(1,0,0,0,0,0,(2.4942,'1')),"reaction affinity":((-1,"A"),0,0,0,0,0,(-2.4942,"A"))}

        return sheet

    def buildDesiredParameters(self):
        '''
        builds up a dictionary of identifiers for those parameters that are chosen to
        be used for the balancing.
        they are of the kind: [(equilibrium constant,reaction)]=keq_1
        '''
        desired_parameters       = {}
        for parameter_type in self.parameter_dict.keys():
            if self.parameter_dict[parameter_type] == True:
                if parameter_type in self.species_parameters:
                    for species in self.species_list:
                        desired_parameters[parameter_type,species] = self.quantity2identifier[parameter_type]+'_'+str(self.species2number[species])
                elif parameter_type in self.reaction_parameters:
                    for reaction in self.reaction_list:
                        desired_parameters[parameter_type,reaction] = self.quantity2identifier[parameter_type]+'_'+str(self.reaction2number[reaction])
                elif parameter_type in self.reaction_species_parameters:
                    for reaction_species in self.model_specific:
                        if parameter_type == reaction_species[0]:
                            desired_parameters[parameter_type,reaction_species[1],reaction_species[2]] = self.quantity2identifier[parameter_type]+'_'+str(self.reaction2number[reaction_species[1]])+'_'+str(self.species2number[reaction_species[2]])

        return desired_parameters

    def collectAvailableValues(self,inconsistent_old_rows):
        '''
        generates the x-vector (values that are provided by the SBtab).
        we have to get the values from the initial SBtab first, afterwards
        the values that have been added by the use of priors and pseudos.
        so we are basically searching two SBtab files and storing their content.
        '''
        self.x_types           = []
        self.parameter2bounds  = {}

        #first: get the column indices of the initial SBtab file, so we can search it
        self.getColumnIndices(self.old_header)
        self.old_rows = self.tidyUpSBtab(inconsistent_old_rows,False,True)

        value_tuples    = []
        row_identifiers = []

        for row in self.old_rows:
            quantity = row[self.qt_column]
            if quantity in self.parameter_dict.keys():
                #x = self.parameter_dict[quantity]    #this check is needed: see, if we know this quantity type
                if len(row) == len(self.old_header):
                    single_tuple = []
                    single_tuple.append(row[self.qt_column])
                    single_tuple.append(row[self.r_column])
                    single_tuple.append(row[self.s_column])
                    single_tuple.append(row[self.m_column])
                    if not row[self.std_column] == '': single_tuple.append(row[self.std_column])
                    elif row[self.qt_column] in self.quantity_type2mean_std.keys(): single_tuple.append('35.0')
                    else: single_tuple.append(str(float(row[self.m_column])*0.5))
                    single_tuple.append(self.makeIdentifier(row))
                    self.x_types.append(single_tuple[5])
                   
                    #check, whether this parameter is really part of the model structure and store boundaries
                    if not row[self.r_column] == '' and not row[self.s_column] == '':
                        row_identifier = (row[self.qt_column],row[self.r_column],row[self.s_column])
                        if self.min_column and self.max_column:
                            self.parameter2bounds[row[self.qt_column],(row[self.r_column],row[self.s_column])] = (row[self.min_column],row[self.max_column])
                        else: self.parameter2bounds[row[self.qt_column],(row[self.r_column],row[self.s_column])] = (None,None)
                    elif not row[self.r_column] == '':
                        row_identifier = (row[self.qt_column],row[self.r_column])
                        if self.min_column and self.max_column:
                            self.parameter2bounds[row[self.qt_column],row[self.r_column]] = (row[self.min_column],row[self.max_column])
                        else: self.parameter2bounds[row[self.qt_column],(row[self.r_column])] = (None,None)
                    else:
                        row_identifier = (row[self.qt_column],row[self.s_column])
                        if self.min_column and self.max_column:
                            self.parameter2bounds[row[self.qt_column],row[self.s_column]] = (row[self.min_column],row[self.max_column])
                        else: self.parameter2bounds[row[self.qt_column],(row[self.s_column])] = (None,None)
                    row_identifiers.append(row_identifier)
                    if row_identifier in self.desired_parameters.keys(): value_tuples.append(single_tuple)
            else:
                raise ParameterBalancingError('There is a quantity type in the SBtab that cannot be interpreted: %s'%(row[self.qt_column]))


        #second: get column indices of the new SBtab file and search it
        self.getColumnIndices(self.new_header)
        for row in self.new_rows:
            if len(row) == len(self.new_header):
                single_tuple = []
                single_tuple.append(row[self.qt_column])
                single_tuple.append(row[self.r_column])
                single_tuple.append(row[self.s_column])
                if not row[self.r_column] == '' and not row[self.s_column] == '':
                    row_identifier = (row[self.qt_column],(row[self.r_column],row[self.s_column]))
                elif not row[self.r_column] == '': row_identifier = (row[self.qt_column],row[self.r_column])
                else: row_identifier = (row[self.qt_column],row[self.s_column])

                if not row_identifier in row_identifiers:
                    if self.min_column and self.max_column:
                        self.parameter2bounds[row_identifier] = (None,None)                
                if row[self.qt_column] in self.prior_list and not self.prior_used: continue
                elif row[self.qt_column] in self.pseudo_list and not self.pseudo_used: continue
                
                if not row_identifier in row_identifiers: # and row[self.m_column] != '':
                    single_tuple.append(row[self.m_column])
                    single_tuple.append(row[self.std_column])
                    single_tuple.append(self.makeIdentifier(row))
                    self.x_types.append(single_tuple[5])
                    if self.min_column and row[self.min_column] != '':
                        single_tuple.append((float(row[self.min_column]),float(row[self.max_column])))
                    elif self.min_column:
                        single_tuple.append((None,None))
                   
                    value_tuples.append(single_tuple)
                    if self.min_column and self.max_column:
                        self.parameter2bounds[row_identifier] = (None,None)
                    #do we have to ask up here if the row_ident is already in parameter2bounds.keys?

        #generate logarithms
        means = []
        stds  = []
        types = []

        for single_tuple in value_tuples:
            means.append(float(single_tuple[3]))
            stds.append(float(single_tuple[4]))
            types.append(single_tuple[0])

        (self.log_means_x,self.log_stds_x) = self.normal2log(means,stds,types)

        return value_tuples

    def makeIdentifier(self,row):
        '''
        generates a specific identifier that works for only this specific row
        '''
        #qidentifier: kG, kV, kM, etc.
        qidentifier = self.quantity2identifier[row[self.qt_column]]
        #qnumber: specific species and/or reaction number in the order of their appearance in SBML model
        if row[self.qt_column] in self.species_parameters: qnumber = str(self.species2number[row[self.s_column]])
        elif row[self.qt_column] in self.reaction_parameters: qnumber = str(self.reaction2number[row[self.r_column]])
        else: qnumber = str(self.species2number[row[self.s_column]])+'_'+str(self.reaction2number[row[self.r_column]])
        identifier = qidentifier+'_'+qnumber

        return identifier

    def getShannons(self):
        '''
        calculate the shannon entropies of the prior and the posterior covariance matrices
        '''
        shannons = []
        #first, shannon entropy of prior matrix self.C_0
        shannon_prior = 0.5 * numpy.log((2*numpy.pi*numpy.exp(2)) * scipy.linalg.det(self.C_0))
        shannons.append(shannon_prior)

        #second, get shannon entropy of posterior matrix
        shannon_posterior = 0.5 * numpy.log((2*numpy.pi*numpy.exp(2)) * scipy.linalg.det(self.C_post_inc))
        shannons.append(shannon_posterior)

        return shannons

    def makeCpostString(self):
        '''
        in order to export the posterior covariance matrix to the user, we make a tsv-string out of it (nicer to read)
        '''
        new_C = '\n'.join(('\t'.join(str(e) for e in row)) for row in self.C_post)

        return new_C

    def getSamplesFromPosterior(self,posterior,posterior_inc,c_matrix_inc,r_matrix,new_SBtab,header,number):
        '''
        give the posterior to this function and receive *number* times a sample from the posterior distribution
        '''
        #first, get the matrix root of the posterior covariance matrix
        C_root = scipy.linalg.matfuncs.sqrtm(c_matrix_inc)
        list_of_SBtab_strings = []

        #do for *number* times:
        for i in range(number):
            #generate random variables of posteriors length
            ksi = numpy.random.normal(0,1,len(posterior_inc))

            #third: get the new posterior for the new SBtab
            new_posterior = numpy.array(posterior)+(numpy.dot(numpy.dot(r_matrix,C_root),ksi))

            #last, but not least: get the new SBtab with the sampled values
            (header,new_SBtab) = self.getNewSBtab(new_sbtab=new_SBtab,posterior_sample=new_posterior,header=header)

            list_of_SBtab_strings.append(new_SBtab)

        return list_of_SBtab_strings 

    def buildDependenceMatrix(self):
        '''
        builds the dependence matrix D from all the submatrices needed for the balancing
        '''
        #the dependence matrix consists of two major parts: the unit matrix on top
        #and the rows below the unit matrix
        D_matrix = []
        self.parameter2row = {}
        
        #first, we build up the unit matrix
        unit_matrix = self.buildUnitMatrix()
        for row in unit_matrix: D_matrix.append(row)

        #second, we check which bottom rows we have to build up
        for pseudo_quantity in self.pseudo_list:
            if self.parameter_dict[pseudo_quantity]:
                #check, which submatrices are needed (are some priors excluded?)
                all_submatrices    = self.sheet[pseudo_quantity]
                needed_submatrices = []
                for i,prior_quantity in enumerate(self.prior_list):
                    if self.parameter_dict[prior_quantity]:
                        needed_submatrices.append((prior_quantity,all_submatrices[i]))
                #check, how many rows the quantity requires and build them via buildBottomRow
                rows = self.buildBottomRow(pseudo_quantity,needed_submatrices)

                for row in rows:
                    D_matrix.append(row)

        matrix = numpy.array(D_matrix)
        return matrix

    def buildUnitMatrix(self):
        '''
        builds up the unit matrix as a first part of the dependence matrix D
        uses only the prior parameters that are chosen by the user
        '''
        unit_rows       = []
        self.type_order = []
        self.bounds     = []
        self.id_order   = {}

        for i,x in enumerate(self.theta_vector):
            row    = [0.0]*len(self.theta_vector)
            row[i] = 1.0
            unit_rows.append(row)
            self.parameter2row[(x[0],x[2])] = row
            self.type_order.append(x[0])
            self.id_order[(x[0],x[2])] = i
            if self.min_column:
                self.bounds.append(self.parameter2bounds[(x[0],x[2])])

        self.matrix_row_counter = len(self.id_order)
        return unit_rows
 
    def buildBottomRow(self,pseudo_quantity,submatrices):
        '''
        builds one of the bottom rows of the dependence matrix D individually.
        of course, this function is not short. but the build up of the specific
        rows of D is extremely complex and shall be able to be generic.
        '''
        use_list            = self.quantity2list[pseudo_quantity]
        sheet               = self.sheet[pseudo_quantity]
        self.remember_links = {}        #remember which parameter belongs to which row (important for generating D_x)
        rows                = []

        row_index = 0

        for i,element in enumerate(use_list):
            row          = [0.0]*len(self.theta_vector)
            column_index = 0
            if self.min_column and self.max_column:
                self.bounds.append(self.parameter2bounds[pseudo_quantity,element])
            for j,matrix_type in enumerate(sheet):
                if self.parameter_dict[self.prior_list[j]]:
                    index_list = self.quantity2list[self.prior_list[j]]
                    #build zero matrix
                    if matrix_type == 0:
                        column_index += (len(index_list))
                    #build unit matrix
                    if matrix_type == 1:
                        for k in range(len(index_list)):
                            row[column_index+row_index] = 1.0
                            column_index += 1
                            break
                        column_index += len(index_list)-1
                    #build specific matrices
                    try:
                        if len(matrix_type)>1:
                            factor = matrix_type[0]
                            matrix = matrix_type[1]
                            #build N, the stoichiometric matrix
                            if matrix == 'A':
                                reactants       = self.reactions_reactants[element][0]
                                stoichiometry_r = self.reactions_reactants[element][1]
                                products        = self.reactions_products[element][0]
                                stoichiometry_p = self.reactions_products[element][1]
                                for species in self.species_list:
                                    for r,reactant in enumerate(reactants):
                                        if species == reactant:
                                            row[column_index] = stoichiometry_r[r]*factor
                                            column_index += 1
                                    for p,product in enumerate(products):
                                        if species == product:
                                            row[column_index] = stoichiometry_p[p]*factor
                                            column_index += 1
                                    if species not in reactants and species not in products:
                                        column_index += 1
                            #build Z, the values for the Michaelis constant coefficients
                            elif matrix == 'Z':
                                #for reaction in self.reaction_list:
                                reactants = self.reactions_reactants[element][0]
                                products  = self.reactions_products[element][0]
                                for michaelis_tuple in self.model_michaelis:
                                    if michaelis_tuple[2] in reactants:
                                        row[column_index] = -1.0
                                        column_index += 1
                                    if michaelis_tuple[2] in products:
                                        row[column_index] = 1.0
                                        column_index += 1
                                    if michaelis_tuple[2] not in reactants and michaelis_tuple[2] not in products:
                                        column_index += 1
                            #build 1, a simple alternative to N
                            elif matrix == '1':
                                for k in range(len(index_list)):
                                    row[column_index+row_index] = 1.0*factor
                                    column_index += 1
                                    break
                                column_index += len(index_list)-1
                    except: pass
            row_index += 1
    
            rows.append(row)
            self.parameter2row[(pseudo_quantity,use_list[i])] = row
            self.type_order.append(pseudo_quantity)
            self.id_order[(pseudo_quantity,element)] = self.matrix_row_counter
            self.matrix_row_counter += 1
            
        return rows

    def buildSpecificDependenceMatrix(self):
        '''
        the matrix D_x is the dependence matrix that holds all rows from D which can
        be also found in the x_vector. Thus, the dependence matrix for the available
        values.
        '''
        rows = []
      
        for i,single_tuple in enumerate(self.x_vector):
            if single_tuple[0] in self.species_parameters:
                rows.append(self.parameter2row[(single_tuple[0],single_tuple[2])])
            elif single_tuple[0] in self.reaction_parameters:
                rows.append(self.parameter2row[(single_tuple[0],single_tuple[1])])
            elif single_tuple[0] in self.reaction_species_parameters:
                rows.append(self.parameter2row[(single_tuple[0],(single_tuple[1],single_tuple[2]))])
            else: print('row identifier not found: ',single_tuple[0],' ',single_tuple[1],' ',single_tuple[2])

        matrix = numpy.array(rows)

        return matrix

    def buildThetaVector(self):
        '''
        generates the theta_vector (default prior means for every parameter in model)
        '''
        self.parameter2row = {}
        theta = []
        means = []
        stds  = []
        types = []
        used_identifiers = []
        self.types_inc   = []
        self.bounds_inc  = []

        for quantity in self.prior_list:
            if self.parameter_dict[quantity]:
                if quantity in self.species_parameters:
                    for species in self.species_list:
                        if not (quantity,species) in used_identifiers:
                            theta.append((quantity,self.prior_values[quantity][0][0],species))
                            means.append(self.prior_values[quantity][0][0])
                            stds.append(self.prior_values[quantity][0][1])
                            types.append(quantity)
                            self.types_inc.append(quantity)
                            if (quantity,species) in self.parameter2bounds.keys():
                                if self.min_column:
                                    self.bounds_inc.append(self.parameter2bounds[(quantity,species)])
                            used_identifiers.append((quantity,species))
                elif quantity in self.reaction_parameters:
                    for reaction in self.reaction_list:
                        if not (quantity,reaction) in used_identifiers:
                            theta.append((quantity,self.prior_values[quantity][0][0],reaction))
                            means.append(self.prior_values[quantity][0][0])
                            stds.append(self.prior_values[quantity][0][1])
                            types.append(quantity)
                            self.types_inc.append(quantity)
                            if (quantity,reaction) in self.parameter2bounds.keys():
                                if self.min_column:
                                    self.bounds_inc.append(self.parameter2bounds[(quantity,reaction)])
                            used_identifiers.append((quantity,reaction))
                elif quantity in self.reaction_species_parameters:
                    for reaction_species in self.model_dict[quantity]:
                        if not (quantity,(reaction_species[1],reaction_species[2])) in used_identifiers:                        
                            theta.append((quantity,self.prior_values[quantity][0][0],(reaction_species[1],reaction_species[2])))
                            means.append(self.prior_values[quantity][0][0])
                            stds.append(self.prior_values[quantity][0][1])
                            types.append(quantity)
                            self.types_inc.append(quantity)
                            if (quantity,(reaction_species[1],reaction_species[2])) in self.parameter2bounds.keys():
                                if self.min_column:
                                    self.bounds_inc.append(self.parameter2bounds[(quantity,(reaction_species[1],reaction_species[2]))])
                            used_identifiers.append((quantity,(reaction_species[1],reaction_species[2])))

        (self.log_means_prior,self.log_stds_prior) = self.normal2log(means,stds,types)

        return theta

    def buildCovarianceMatrices(self):
        '''
        generate covariance matrix for measured values x (stds) from SBtab
        '''
        #first, generate prior covariance matrix C_0
        C_0_rows = []
        for i,theta in enumerate(self.theta_vector):
            row = [0.0]*len(self.theta_vector)
            try: row[i] = numpy.square(float(self.getDefaultStd(self.theta_vector[i][0])))
            except: raise ParameterBalancingError("The prior covariance matrix cannot be built satisfactory.")
            C_0_rows.append(row)

        C_0 = numpy.array(C_0_rows)

        #second, generate covariance matrix according to the input values in the x-vector
        C_x_rows = []
        for i,x_entry in enumerate(self.x_vector):
            row = [0.0]*len(self.x_vector)
            if self.log_stds_x[i] == 0.0 or self.log_stds_x[i] == 'NaN' or numpy.square(self.log_stds_x[i]) == 0.0 or self.log_stds_x[i]<0.00001:
                row[i] = numpy.square(0.0001)
            else:
                row[i] = numpy.square(self.log_stds_x[i])
            ####
            #if x_entry[0] == 'equilibrium constant':
            #    row[i] = numpy.square(0.0001)
            ####
            C_x_rows.append(row)

        C_x = numpy.array(C_x_rows)
        
        return C_0, C_x


    def getDefaultStd(self,quantity):
        '''
        returns the default standard deviation
        '''
        return self.log_stds[name2index[quantity]]


    def getDefaultMean(self,quantity):
        '''
        returns the default mean value of quantity
        '''
        return self.means[name2index[quantity]]

    def tryStoringConditions(self):
        '''
        if the SBtab file holds information on the measurement conditions (temp or pH), store them
        '''
        self.stored_temp = []
        self.stored_pH   = []

        checklength = len(self.header)
        for i,row in enumerate(self.rows):
            if not len(row) != checklength:
                if self.temp_column:
                    stored_temp = []
                    stored_temp.append(row[self.qt_column])
                    stored_temp.append(row[self.r_column])
                    stored_temp.append(row[self.s_column])
                    stored_temp.append(row[self.temp_column])
                    self.stored_temp.append(stored_temp)
                if self.pH_column:
                    stored_pH = []
                    stored_pH.append(row[self.qt_column])
                    stored_pH.append(row[self.r_column])
                    stored_pH.append(row[self.s_column])
                    stored_pH.append(row[self.pH_column])
                    self.stored_pH.append(stored_pH)

    def calculatePosteriori(self):
        '''
        calculates the posteriori values
        '''
        #posterior covariance
        try:
            self.C_0_inv = numpy.linalg.inv(self.C_0)
        except:
            raise ParameterBalancingError("C_0 is not invertible\n"+str(self.C_0.shape)+"\n"+str(self.C_0)+"\n"+str(e))
        try:
            self.C_x_inv = numpy.linalg.inv(self.C_x)
        except:
            raise ParameterBalancingError("C_x is not invertible\n"+str(self.C_x.shape)+"\n"+str(self.C_x)+"\n"+str(e))

        self.C_post_inc = numpy.linalg.inv(self.C_0_inv+numpy.dot(numpy.dot(self.D_x.transpose(),self.C_x_inv),self.D_x))
        self.C_post = numpy.dot((numpy.dot(self.D,self.C_post_inc)),self.D.transpose())

        #posterior vector
        self.mean_post_inc = numpy.dot(self.C_post_inc,(numpy.dot(numpy.dot(self.D_x.transpose(),self.C_x_inv),self.log_means_x)+numpy.dot(self.C_0_inv,self.log_means_prior)))

        self.means_log_inc = copy.deepcopy(self.mean_post_inc)
        self.stds_log_inc  = self.extractPostCovarianceInc()
        
        self.log_mean_post = numpy.dot(self.D,self.mean_post_inc)
                
    def extractPostCovariance(self):
        '''
        extract the stds from the posterior diagonal covariance matrix C_post
        '''
        stds = []
        for i,row in enumerate(self.C_post):
            #if float(row[i])<0.: stds.append(row[i])
            #else: stds.append(numpy.sqrt(row[i]))
            stds.append(numpy.sqrt(row[i]))
        return stds
    
    def extractPostCovarianceInc(self):
        '''
        extract the stds from the posterior diagonal covariance matrix C_post
        '''
        stds = []
        for i,row in enumerate(self.C_post_inc):
            #if name2index[self.types_inc[i-1]] in thermodynamic_indices: stds.append(row[i])
            #else: stds.append(numpy.sqrt(row[i]))
            stds.append(numpy.sqrt(row[i]))
        return stds

    def getNewSBtab(self,new_sbtab=False,posterior_sample=False,header=False):
        '''
        generates new SBtab
        '''
        if new_sbtab:
            self.new_rows = []
            for row in new_sbtab: self.new_rows.append(row)
            means           = posterior_sample
            try: self.new_header = header.split('\t')
            except: self.new_header = header
        else:
            means = self.mean_post

        if not self.medianed: means = self.old_means
        else: means = self.mean_post_opt

        finished_rows = [self.new_header]
        first         = True
        self.hilo     = []

        # Add biomass warning if required
        self.check4biomass()

        for i,row in enumerate(self.new_rows):
            if len(row) == len(self.new_header):
                if row[0] == 'QuantityType': continue
                #first: identify the row
                if row[0] in self.species_parameters: row_identifier = (row[0],row[2])
                elif row[0] in self.reaction_parameters: row_identifier = (row[0],row[1])
                else: row_identifier = (row[0],(row[1],row[2]))
                try: row_number = self.id_order[row_identifier]
                except: continue
                #second: fill the row with the balanced values
                row[5] = str(format(float(means[row_number]),'.4f'))
                row[6] = str(format(float(self.stds_post[row_number]),'.4f'))
                if row[0] in self.quantity_type2mean_std.keys():
                    row[3] = str(format(float(means[row_number]),'.4f'))
                else:
                    (log_mean,log_std) = self.normal2log([row[5]],[row[6]],row[0])
                    row[3] = str(format(float(numpy.exp(log_mean[0])),'.4f'))
                if not row[0] in thermodynamics:
                    row[self.logm_column] = str(format(float(self.log_mean_post[row_number]),'.4f'))
                    row[self.logs_column] = str(format(float(self.log_stds_post[row_number]),'.4f'))
                finished_rows.append(row)
                first = self.checkHighAndLow(row,first)
        if not first: self.log += '\n'
        if self.hilo != []:
            self.hilo = sorted(self.hilo)
            for entry in self.hilo:
                self.log += entry

        return self.new_header,finished_rows

    def checkHighAndLow(self,row,first):
        '''
        this function checks whether the given parameter (in form of its corresponding SBtab row)
        has a significantly high or low value.
        '''
        if self.pmin[row[0]] == None and self.pmax[row[0]] == None: return first
        elif float(row[5]) < self.pmin[row[0]]:
            if first:
                self.log += '\n### Warnings about unusually high or low values ###\n'
                first = False
            #self.log += 'The value for the %s of %s,%s lies under the given lower bound: %s. Please check the accuracy and refer to the FAQ for help.\n'%(row[0],row[1],row[2],row[4])
            self.hilo.append('The value for the %s of %s,%s lies under the given lower bound: %s. Please check the accuracy and refer to the FAQ for help.\n'%(row[0],row[1],row[2],row[5]))
        elif float(row[5]) > self.pmax[row[0]]:
            if first:
                self.log += '### Warnings about unusually high or low values ###\n'
                first = False
            #self.log += 'The value for the %s of %s,%s lies over the given upper bound: %s. Please check the accuracy and refer to the FAQ for help.\n'%(row[0],row[1],row[2],row[4])
            self.hilo.append('The value for the %s of %s,%s lies over the given upper bound: %s. Please check the accuracy and refer to the FAQ for help.\n'%(row[0],row[1],row[2],row[5]))

        return first

    def mimicMean(self,badSBtab):
        '''
        since Wolf wants only the MEAN to be displayed instead of the VALUE, but the kineticizer really
        needs the VALUE, this function renames it
        '''
        header       = badSBtab.split('\n')[0]
        new_header   = ''
        splitter     = header.split('\t')

        for i,column in enumerate(splitter):
            if self.medianed:
                if column == 'Mean':
                    new_header += 'Value\t'
                    v_column = i
                elif column == 'Std':
                    new_header += str(column)+'\t'
                    s_column = i
                elif column == 'Mode':
                    new_header += str(column)+'\t'
                    med_column = i               
                else:
                    new_header += str(column)+'\t'
            else:
                if column == 'Mode':
                    new_header += 'Value\t'
                    v_column = i
                elif column == 'Std':
                    new_header += str(column)+'\t'
                    s_column = i
                elif column == 'Mean':
                    new_header += str(column)+'\t'
                    med_column = i
                else:
                    new_header += str(column)+'\t'

        x_header  = str(new_header[:-1]) + '\n'
        goodSBtab = x_header

        for row in badSBtab.split('\n')[1:]:
            no_tabs = row.split('\t')
            if len(no_tabs) == len(splitter) and not no_tabs[0].startswith('standard') and not no_tabs[0].startswith('chemical') and not no_tabs[0].startswith('reaction'):
                x_row   = ''
                new_row = ''
                for value in no_tabs:
                    x_row  += str(value)+'\t'
                    new_row = str(x_row[:-1])
                goodSBtab  += new_row+'\n'
            elif no_tabs[0].startswith('standard') or no_tabs[0].startswith('standard') or no_tabs[0].startswith('chemical') or no_tabs[0].startswith('reaction'):
                x_row   = ''
                new_row = ''
                for value in no_tabs:
                    x_row += str(value)+'\t'
                new_row    = str(x_row[:-1])
                goodSBtab += new_row+'\n'

        return goodSBtab

    def checkForBlanks(self,oldSBtab):
        '''
        if the SBtab file shall be inserted into the SBML file, every (!) value field has to hold a value.
        this is checked by this function
        '''
        tabbed_rows = oldSBtab.split('\n')[1:]
        headertab   = oldSBtab.split('\n')[0]
        newrows     = headertab+'\n'
        header      = headertab.split('\t')
		
        for rowtab in tabbed_rows:
            row = rowtab.split('\t')
            if len(row) == len(header):
                if row[5] == '' or row[5] == 'nan':
                    row[5] = '1'
                    if row[1] != '' and row[2] == '':
                        print('Could not find a value for %s of reaction %s. Value set to 1.'%(row[0],row[1]))
                    elif row[1] == '' and row[2] != '':
                        print('Could not find a value for %s of species %s. Value set to 1.'%(row[0],row[2]))
                    elif row[1] != '' and row[2] != '':
                        print('Could not find a value for %s of reaction %s and species %s. Value set to 1.'%(row[0],row[1],row[2]))
                    else:
                        print('Could not find a value for a specific kinetic parameter. Value set to 1.')
            newrow = ''
            for value in row:
                newrow += value
                newrow += '\t'
            newrows += newrow[:-1]+'\n'

        return newrows



    def expand_R_x_theta(self):
        '''
        if we get the informations of pH and temperature, we have to adjust the dimension of R_x_theta
        by simply adding some 0-columns
        '''
        for row in self.matrix_array:
            for i in range(len(self.species_list)*2):
                row.insert(len(self.species_list)-len(self.enzymeStorage)+i,0.0)

        self.R_x_theta = numpy.array(self.matrix_array) #scipy.sparse.csc_matrix(self.matrix_array) #numpy.array(self.matrix_array)


    def get_R_0(self):
        '''
        generate R-matrix for the values that are provided; it doesnt need to be complete
        '''
        rows               = []
        rows_mu_t          = []
        rows_mu_pH         = []
        rows_eq_t          = []
        rows_eq_pH         = []
        rows_mu_t_name     = []
        rows_mu_pH_name    = []
        rows_eq_t_name     = []
        rows_eq_pH_name    = []        
        names              = []
        x_names            = []
        self.parameter2row = {}
        #self.delta_temp    = []
        #self.delta_pH      = []
        #self.delta_eq_temp = []
        #self.delta_eq_pH   = []

        for entry in self.x_vector:
            if entry[0] in self.x_all:
                x_names.append(entry[0])

        R_x_theta_dict = {}

        for i,name in enumerate(self.x_all):
            R_x_theta_dict[name] = self.matrix_array[i]

        for i,parameter_name in enumerate(x_names):
            self.parameter2row[parameter_name] = R_x_theta_dict[parameter_name]

        for i,available_parameter in enumerate(self.x_vector):
            for i,parameter in enumerate(available_parameter):
                if available_parameter[i] == 'nan':
                    available_parameter[i] = ''
                    
            if available_parameter[0] in self.x_all:
                row  = copy.deepcopy(self.parameter2row[available_parameter[0]])
                name = available_parameter[0]
                #if not name in names:

                # for regression: check whether there are temperatures and pH values
                # if so, calculate the delta_temp and delta_pH, further calculate the
                # matrices for the additional columns (mu_t and mu_pH)
                if self.temp_column and str(available_parameter[0]).startswith('01_k_G'):
                    row_identifier = []
                    row_identifier.extend([available_parameter[5],available_parameter[1],available_parameter[2]])
                    for single_temp_store in self.stored_temp:
                        temp_identifier = []
                        temp_identifier.extend([single_temp_store[0],single_temp_store[1],single_temp_store[2]])
                        if row_identifier == temp_identifier:
                            delta_temp = self.get_delta(single_temp_store[3],'temp')
                            row_copy = copy.deepcopy(row[:len(self.species_list)])
                            row_mu_t = []
                            for number in row_copy:
                                row_mu_t.append(float(number)*float(delta_temp))
                            rows_mu_t.append(row_mu_t)
                            rows_mu_t_name.append(available_parameter[0])
                            break

                if self.temp_column and str(available_parameter[0]).startswith('08_k_eq'):
                    #rows_eq_t      = []
                    row_identifier = []
                    row_identifier.extend([available_parameter[5],available_parameter[1],available_parameter[2]])
                    for single_temp_store in self.stored_temp:
                        if single_temp_store[0].startswith('equilibrium'):
                            temp_identifier = []
                            temp_identifier.extend([single_temp_store[0],single_temp_store[1],single_temp_store[2]])
                            if row_identifier == temp_identifier:
                                delta_eq_temp = self.get_delta(single_temp_store[3],'temp')
                                row_copy = copy.deepcopy(row[:len(self.species_list)])
                                row_eq_t = []
                                for number in row_copy:
                                    row_eq_t.append(float(number)*float(delta_eq_temp))
                                rows_eq_t.append(row_eq_t)
                                rows_eq_t_name.append(available_parameter[0])
                                break

                if self.pH_column and str(available_parameter[0]).startswith('01_k_G'):
                    row_identifier = []
                    row_identifier.extend([available_parameter[5],available_parameter[1],available_parameter[2]])
                    for single_pH_store in self.stored_pH:
                        if single_pH_store[0].startswith('standard'):
                            pH_identifier = []
                            pH_identifier.extend([single_pH_store[0],single_pH_store[1],single_pH_store[2]])
                            if row_identifier == pH_identifier:
                                delta_pH = self.get_delta(single_pH_store[3],'pH')

                                row_copy = copy.deepcopy(row[:len(self.species_list)])
                                row_mu_pH = []
                                for number in row_copy:
                                    row_mu_pH.append(float(number)*float(delta_pH))
                                rows_mu_pH.append(row_mu_pH)
                                rows_mu_pH_name.append(available_parameter[0])
                                break


                if self.pH_column and str(available_parameter[0]).startswith('08_k_eq'):
                    #rows_eq_pH     = []
                    row_identifier = []
                    row_identifier.extend([available_parameter[5],available_parameter[1],available_parameter[2]])
                    for single_pH_store in self.stored_pH:
                        if single_pH_store[0].startswith('equilibrium'):
                            pH_identifier = []
                            pH_identifier.extend([single_pH_store[0],single_pH_store[1],single_pH_store[2]])
                            if row_identifier == pH_identifier:
                                delta_eq_pH = self.get_delta(single_pH_store[3],'pH')
                                row_copy  = copy.deepcopy(row[:len(self.species_list)])
                                row_eq_pH = []
                                for number in row_copy:
                                    row_eq_pH.append(float(number)*float(delta_eq_pH))
                                rows_eq_pH.append(row_eq_pH)
                                rows_eq_pH_name.append(available_parameter[0])
                                break

            rows.append(row)
            names.append(name)

        
        #insert the new columns for mu_temp and mu_pH (only in rows k_G and k_eq)
        if rows_mu_t and rows_mu_pH:
            eq_counter = 0
            for i,row in enumerate(rows):
                if names[i].startswith('01_k'):
                    if rows_mu_t_name.count(rows_mu_t_name[eq_counter])>1 or sum([1 for x in rows_mu_t if len(x)>1])>0:
                        for j in range(len(rows_mu_t[eq_counter])):
                            row[len(self.species_list)-len(self.enzymeStorage)+j] = rows_mu_t[eq_counter][j]
                    else:
                        for j in range(len(rows_mu_t[eq_counter])):
                            row[len(self.species_list)-len(self.enzymeStorage)+j] = 0.0
                    eq_counter += 1
                    if eq_counter == len(rows_eq_t_name):
                        break
                            
            eq_counter = 0
            for i,row in enumerate(rows):
                if names[i].startswith('01_k'):
                    if rows_mu_pH_name.count(rows_mu_pH_name[eq_counter])>1 or sum([1 for x in rows_mu_pH if len(x)>1])>0:
                        for j in range(len(rows_mu_pH[eq_counter])):
                            row[len(self.species_list)*2-len(self.enzymeStorage)*2+j] = rows_mu_pH[eq_counter][j]
                    else:
                        for j in range(len(rows_mu_pH[eq_counter])):
                            row[len(self.species_list)*2-len(self.enzymeStorage)*2+j] = 0.0
                    eq_counter += 1
                    if eq_counter == len(rows_eq_t_name):
                        break


        if rows_eq_t and rows_eq_pH:
            eq_counter = 0
            for i,row in enumerate(rows):
                if names[i].startswith('08_'):
                    if rows_eq_t_name.count(rows_eq_t_name[eq_counter])>1 or sum([1 for x in rows_eq_t if len(x)>1])>0:
                        for j in range(len(rows_eq_t[eq_counter])):
                            row[len(self.species_list)-len(self.enzymeStorage)+j] = rows_eq_t[eq_counter][j]
                    else:
                        for j in range(len(rows_eq_t[eq_counter])):
                            row[len(self.species_list)-len(self.enzymeStorage)+j] = 0.0
                    if rows_eq_pH_name.count(rows_eq_pH_name[eq_counter])>1 or sum([1 for x in rows_eq_pH if len(x)>1])>0:
                        for j in range(len(rows_eq_pH[eq_counter])):
                            row[len(self.species_list)*2-len(self.enzymeStorage)*2+j] = rows_eq_pH[eq_counter][j]
                    else:
                        for j in range(len(rows_eq_pH[eq_counter])):
                            row[len(self.species_list)*2-len(self.enzymeStorage)*2+j] = 0.0
                    eq_counter += 1
                    if eq_counter == len(rows_eq_t_name):
                        break

        R_0 = numpy.array(rows)

        return R_0

    def get_delta(self,stored_value,mode):
        '''
        calculates the delta-value for temperature and pH (measured value minus wanted value)
        '''
        if mode == 'temp':
            if stored_value == 'nan' or stored_value == '':
                stored_value = self.temperature
            delta_value = float(self.temperature)-float(stored_value)
            #rint delta_value,'t'
        elif mode == 'pH':
            if stored_value == 'nan' or stored_value == '':
                stored_value = self.pH
            delta_value = float(self.pH)-float(stored_value)
            #rint delta_value

        return delta_value
