#!/usr/bin/env python
try: from . import SBtab
except: import SBtab
import numpy
import scipy.linalg
import copy
import time
import datetime
import os
import sys

# common declarations for the webversion
name2index = {'standard chemical potential': 0,
              'catalytic rate constant geometric mean': 1,
              'Michaelis constant': 2,
              'activation constant': 3,
              'inhibitory constant': 4,
              'concentration of enzyme': 5,
              'concentration': 6,
              'equilibrium constant': 7,
              'substrate catalytic rate constant': 8,
              'product catalytic rate constant': 9,
              'forward maximal velocity': 10,
              'reverse maximal velocity': 11,
              'chemical potential': 12,
              'reaction affinity': 13}
header_names = ['!QuantityType', '!Reaction:SBML:reaction:id',
                '!Compound:SBML:species:id', '!Mode', '!Unit',
                '!UnconstrainedGeometricMean', '!UnconstrainedGeometricStd',
                '!UnconstrainedMean', '!UnconstrainedStd']
inhibitory_sbos = [20, 206, 207, 536, 537]
activation_sbos = [13, 21, 459, 461, 462]


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
    def __init__(self, sbml_model, req=True):
        '''
        initialise pb class
        '''
        self.model = sbml_model

        # initialise log file
        self.log = 'Parameter balancing log file; %s\n\n'%(time.asctime())

        # rudimentary validity check and model initialisation
        if req:
            try: self.model.getListOfSpecies()
            except:
                raise ParameterBalancingError('You have not used a valid SBML model.')
            self.gain_model_information()

    def check_biomass(self):
        '''
        if the model appears to have a biomass function, we add a warning to
        the log file; biomass reactions should not have a reversible convenience
        kinetic and should be treated differently.
        '''
        for species in self.model.getListOfSpecies():
            if 'biomass' in species.getName().lower():
                self.log += '### Warning about apparent biomass reaction: '\
                            '### \nIt appears that the model has a biomass '\
                            'reaction (SBML ID %s). The convenience kinetics '\
                            'that are entered for the model reactions are not '\
                            'appropriate for biomass reactions. We suggest '\
                            'you replace the kinetics with something more '\
                            'appropriate like the one proposed by Hofmeyr et '\
                            'al. (2013). \n\n' % (species.getId())

        for reaction in self.model.getListOfReactions():
            if 'biomass' in reaction.getName().lower():
                self.log += '### Warning about apparent biomass reaction: '\
                            '### \nIt appears that the model has a biomass '\
                            'reaction (SBML ID %s). The convenience kinetics '\
                            'that are entered for the model reactions are not '\
                            'appropriate for biomass reactions. We suggest '\
                            'you replace the kinetics with something more '\
                            'appropriate like the one proposed by Hofmeyr et '\
                            'al. (2013). \n\n' % (reaction.getId())
                continue

            if reaction.getNumReactants() > 4:
                self.log += '### Warning about apparent biomass reaction: '\
                            '### \nThere is a reaction with an usually high '\
                            'amount of reactants (SBML ID %s). It may be a '\
                            'biomass reaction. The convenience kinetics that '\
                            'are entered for the model reactions are not '\
                            'appropriate for these reactions. We suggest you '\
                            'replace the kinetics with something more '\
                            'appropriate like the one proposed by Hofmeyr et '\
                            'al. (2013). \n\n''' % (reaction.getId())

    def check_enzyme_species(self):
        '''
        sometimes modellers provide enzymes as SBML species, which does not
        really comply to modelling standards; parameter balancing cannot
        correctly assign these enzymes as reaction modifiers under this
        circumstance; thus, we ignore them and add the information to the log.
        for this, we search for modifiers which do not appear in any reaction
        as reactant or product.
        '''
        involved_species = []
        modifiers = []
        
        for reaction in self.model.getListOfReactions():
            for species in reaction.getListOfReactants():
                involved_species.append(species.getId())
            for species in reaction.getListOfProducts():
                involved_species.append(species.getId())
            for species in reaction.getListOfModifiers():
                modifiers.append(species.getSpecies())

        for m_id in modifiers:
            if not m_id in involved_species:
                self.log += ('The modifier %s is not involved in any reaction'\
                             ' as either reactant or product. Please check if'\
                             ' this modifier is required to be an SBML specie'\
                             's. E.g., enzymes should not be SMBL '\
                             'species.' % (m_id))

    def gain_model_information(self):
        '''
        gain all initial information about the SBML model:
        which species and reactions are contained? which modifiers?
        '''
        self.get_parameter_information()

        self.species_list = []
        self.reaction_list = []
        self.species2number = {}
        self.reaction2number = {}

        # collect species and reactions from the SBML model
        for i, species in enumerate(self.model.getListOfSpecies()):
            self.species_list.append(species.getId())
            self.species2number[species.getId()] = i + 1

        for i, reaction in enumerate(self.model.getListOfReactions()):
            self.reaction_list.append(reaction.getId())
            self.reaction2number[reaction.getId()] = i + 1

        (self.model_michaelis,
         self.model_inhibition,
         self.model_activation) = self.check_if_required()
        
        self.model_dict = {'Michaelis constant': self.model_michaelis,
                           'inhibitory constant': self.model_inhibition,
                           'activation constant': self.model_activation}
        
        self.build_reaction_specifics()

        self.quantity2list = {'standard chemical potential': self.species_list,
                              'chemical potential': self.species_list,
                              'concentration': self.species_list,
                              'catalytic rate constant geometric mean': self.reaction_list,
                              'concentration of enzyme': self.reaction_list,
                              'product catalytic rate constant': self.reaction_list,
                              'substrate catalytic rate constant': self.reaction_list,
                              'equilibrium constant': self.reaction_list,
                              'forward maximal velocity': self.reaction_list,
                              'reverse maximal velocity': self.reaction_list,
                              'reaction affinity': self.reaction_list,
                              'Michaelis constant': self.model_michaelis,
                              'inhibitory constant': self.model_inhibition,
                              'activation constant': self.model_activation}

        # add some information for the log file
        self.log += '### Model information ###\n'
        self.log += 'The SBML model has a total of %s Reaction/s and %s '\
                    'Species.\n' % (self.model.getNumReactions(),
                                    self.model.getNumSpecies())

    def get_parameter_information(self):
        '''
        read a table file from the resources directory holding numerous
        informations on the handling of the parameters, the parameter details,
        and how to build the dependency matrix for the different parameter
        types
        '''
        p = os.path.dirname(os.path.abspath(__file__)) + '/files/default_'\
            'files/pb_prior.tsv'
        try: pf = open(p, 'r')
        except:
            print('The prior file (/files/default_files/pb_prior.tsv) coul'\
                  'd not be found. I quit.')
            sys.exit()

        p = pf.read()
        sbtab_prior = SBtab.SBtabTable(p, 'pb_prior.tsv')
 
        self.quantity2identifier = {}
        self.quantity_type2unit = {}
        self.quantity_type2mean_std = {}
        self.quantity_type2median_std = {}
        self.species_parameters = []
        self.reaction_parameters = []
        self.reaction_species_parameters = []
        self.thermodynamics = []
        self.matrix_info = {}
        self.data_std = {}
        prior_list = []
        pseudo_list = []
        
        for row in sbtab_prior.value_rows:
            if not row[sbtab_prior.columns_dict['!QuantityType']] in name2index.keys():
                continue
            if row[sbtab_prior.columns_dict['!UseAsPriorInformation']] == '0':
                continue
            if len(row) == len(sbtab_prior.columns):
                self.quantity2identifier[row[sbtab_prior.columns_dict['!QuantityType']]] = row[sbtab_prior.columns_dict['!Abbreviation']]
                self.quantity_type2unit[row[sbtab_prior.columns_dict['!QuantityType']]] = row[sbtab_prior.columns_dict['!Unit']]

                if row[sbtab_prior.columns_dict['!MathematicalType']] == 'Additive':
                    self.quantity_type2mean_std[row[sbtab_prior.columns_dict['!QuantityType']]] = [float(row[sbtab_prior.columns_dict['!PriorMedian']]),
                                                                                                   float(row[sbtab_prior.columns_dict['!PriorStd']])]
                elif row[sbtab_prior.columns_dict['!MathematicalType']] == 'Multiplicative':
                    self.quantity_type2median_std[row[sbtab_prior.columns_dict['!QuantityType']]] = [float(row[sbtab_prior.columns_dict['!PriorMedian']]),
                                                                                                     float(row[sbtab_prior.columns_dict['!PriorGeometricStd']])]

                if row[sbtab_prior.columns_dict['!Dependence']] == 'Basic':
                    prior_list.append(row[sbtab_prior.columns_dict['!QuantityType']])
                elif row[sbtab_prior.columns_dict['!Dependence']] == 'Derived':
                    pseudo_list.append(row[sbtab_prior.columns_dict['!QuantityType']])

                if row[sbtab_prior.columns_dict['!BiologicalElement']] == 'Species':
                    self.species_parameters.append(row[sbtab_prior.columns_dict['!QuantityType']])
                elif row[sbtab_prior.columns_dict['!BiologicalElement']] == 'Reaction':
                    self.reaction_parameters.append(row[sbtab_prior.columns_dict['!QuantityType']])
                elif row[sbtab_prior.columns_dict['!BiologicalElement']] == 'Reaction/Species':
                    self.reaction_species_parameters.append(row[sbtab_prior.columns_dict['!QuantityType']])

                if row[sbtab_prior.columns_dict['!Unit']] == 'kJ/mol':
                    self.thermodynamics.append(row[sbtab_prior.columns_dict['!QuantityType']])

                self.matrix_info[row[sbtab_prior.columns_dict['!QuantityType']]] = row[sbtab_prior.columns_dict['!MatrixInfo']]
                self.data_std[row[sbtab_prior.columns_dict['!QuantityType']]] = row[sbtab_prior.columns_dict['!DataStd']]

        self.prior_list = self.sort_list(prior_list)
        self.pseudo_list = self.sort_list(pseudo_list)

    def check_if_required(self):
        '''
        not every reaction/species-combo needs an activation, inhibitory, or
        Michaelis constant. Check which combos DO need them.
        '''
        Michaelis = []
        inhibitors = []
        activators = []
        
        for reaction in self.model.getListOfReactions():
            # first: inhibitory and activation constants
            # (only identifiable by SBO terms)
            for modifier in reaction.getListOfModifiers():
                if modifier.getSBOTerm() in inhibitory_sbos:
                    inhibitors.append(['inhibitory constant', reaction.getId(),
                                       modifier.getSpecies()])
                elif modifier.getSBOTerm() in activation_sbos:
                    activators.append(['activation constant', reaction.getId(),
                                       modifier.getSpecies()])

            #second: Michaelis constants
            for reactant in reaction.getListOfReactants():
                Michaelis.append(['Michaelis constant', reaction.getId(),
                                  reactant.getSpecies()])
            for product in reaction.getListOfProducts():
                Michaelis.append(['Michaelis constant', reaction.getId(),
                                  product.getSpecies()])

        self.model_specific = Michaelis+inhibitors+activators
        return Michaelis, inhibitors, activators

    def build_reaction_specifics(self):
        '''
        generate dictionaries linking the reactions with its reactants,
        products, and stoichiometric coefficients
        '''
        self.reactions_reactants = {}
        self.reactions_products = {}
        
        for i, reaction in enumerate(self.model.getListOfReactions()):
            reactants = []
            stoich = []
            for reactant in reaction.getListOfReactants():
                reactants.append(reactant.getSpecies())
                this_stoich = reactant.getStoichiometry()
                stoich.append(this_stoich*(-1.0))
            self.reactions_reactants[reaction.getId()] = (reactants, stoich)

            products = []
            stoich = []
            for product in reaction.getListOfProducts():
                products.append(product.getSpecies())
                this_stoich = product.getStoichiometry()
                stoich.append(this_stoich)
            self.reactions_products[reaction.getId()] = (products, stoich)

    def make_empty_sbtab(self, pmin, pmax, parameter_dict):
        '''
        if no SBtab is given, create an empty SBtab for the model using
        function make_sbtab and some default parameters
        '''
        sbtab_string = '!!SBtab TableType="Quantity" Version="0.1" Level="1.0" '\
                       'TableName="EmptyParameterFile"\n'\
                       '!QuantityType\t!Reaction:SBML:reaction:id\t'\
                       '!Compound:SBML:species:id\t!Mean\t!Std\t!Unit\n'
        empty_sbtab = SBtab.SBtabTable(sbtab_string, 'empty.csv')
        value_rows = self.make_sbtab(empty_sbtab, 'empty.csv',
                                     'All organisms', 43, pmin, pmax,
                                     parameter_dict)
        return value_rows

    
    def make_sbtab(self, sbtab, file_name, organism, volume,
                   pmin, pmax, parameter_dict):
        '''
        makes an insertable SBtab-file out of a possibly erraneous SBtab-File
        '''        
        # set warning message for numerical problems in numpy
        numpy.seterrcall(self.print_warning)
        numpy.seterr(all='call')
        self.warned = True
        
        self.sbtab = sbtab
        self.organism = organism
        self.new_header = header_names
        self.new_rows = []
        self.pmin = pmin
        self.pmax = pmax
        self.parameter_dict = parameter_dict
        # the possibly messy user file needs to be tidied before computation
        self.rows = self.tidy_up_sbtab(False)
        
        # build a list of all parameters provided by the user
        self.available_parameters = []
        for i, row in enumerate(self.rows):
            if len(row) == len(sbtab.columns):

                if row[self.sbtab.columns_dict['!QuantityType']] in self.species_parameters:
                    self.available_parameters.append([row[self.sbtab.columns_dict['!QuantityType']],
                                                      row[self.sbtab.columns_dict['!Compound:SBML:species:id']]])
                elif row[self.sbtab.columns_dict['!QuantityType']] in self.reaction_parameters:
                    self.available_parameters.append([row[self.sbtab.columns_dict['!QuantityType']],
                                                      row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']]])
                elif row[self.sbtab.columns_dict['!QuantityType']] in self.reaction_species_parameters:
                    self.available_parameters.append([row[self.sbtab.columns_dict['!QuantityType']],
                                                      row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']],
                                                      row[self.sbtab.columns_dict['!Compound:SBML:species:id']]])

        # build all required model parameters with respect to the
        # provided parameters: either collect or create the parameter
        for species in self.species_list:
            for quantity in self.species_parameters:
                if [quantity, species] in self.available_parameters:
                    self.new_rows.append(self.existing_row(species, quantity))
                else:
                    self.new_rows.append(self.new_row(species, quantity))

        for reaction in self.reaction_list:
            for quantity in self.reaction_parameters:
                if [quantity, reaction] in self.available_parameters:
                    self.new_rows.append(self.existing_row(reaction, quantity))
                else:
                    self.new_rows.append(self.new_row(reaction, quantity))

        tuple_list = [self.model_michaelis, self.model_activation,
                      self.model_inhibition]
        for i, t_list in enumerate(tuple_list):
            for reaction_species in t_list:
                if reaction_species in self.available_parameters:
                    self.new_rows.append(self.existing_row(reaction_species,
                                                           self.reaction_species_parameters[i]))
                else:
                    self.new_rows.append(self.new_row(reaction_species,
                                                      self.reaction_species_parameters[i]))

        self.log += 'You have used %s values to describe %s unique kinetic '\
                    'parameters as input for the balancing.\n' % (len(self.rows),
                                                                  len(self.new_rows))
        nr = self.model.getNumReactions()
        self.log += 'The parameter balancing will output %s standard chemical '\
                    'potentials, %s catalytic rate constants geometric mean, '\
                    '%s activation constants, %s inhibitory constants, %s con'\
                    'centrations, %s concentrations of enzyme, %s equilibrium'\
                    'constants, %s substrate catalytic rate constants, %s prod'\
                    'uct catalytic rate constants, %s forward maximal velociti'\
                    'es, %s reverse maximal velocities, %s chemical potentials'\
                    ', and %s reaction affinities.\n'%(self.model.getNumSpecies(),
                                                       nr, nr, nr,
                                                       self.model.getNumSpecies(),
                                                       nr, nr, nr, nr, nr, nr,
                                                       self.model.getNumSpecies(),
                                                       nr)

        # create new SBtab file from the collected information
        sbtab_string = '!!SBtab TableType="Quantity" Version="0.1" Level="1.0" '\
                       'TableName="%s"\n' % (file_name) + '\t'.join(self.new_header) + '\n'

        for row in self.new_rows:
            sbtab_string += '\t'.join(row) + '\n'
        new_sbtab = SBtab.SBtabTable(sbtab_string, file_name)

        return new_sbtab
                    
    def tidy_up_sbtab(self, ig_bounds):
        '''
        remove a lot of stuff from SBtab that is either cramping the functionality or
        delete entries that must not be used in parameter balancing
        '''
        consistent_rows = []
        log_header = False

        try: mean_column = self.sbtab.columns_dict['!Mean']
        except:
            mean_column = self.sbtab.columns_dict['!UnconstrainedGeometricMean']
        try: std_column = self.sbtab.columns_dict['!Std']
        except:
            std_column = self.sbtab.columns_dict['!UnconstrainedGeometricStd']
        
        for row in self.sbtab.value_rows:
            if len(row) == len(self.sbtab.value_rows[0]):
                # irregular entries are set to empty string
                for i, value in enumerate(row):
                    if value == 'nan' or value == '-' or value == 'NaN':
                        row[i] = ''
                    elif value == 'None' or value == None:
                        row[i] = ''
                
                # if the user has specified an organism then remove all others
                if self.organism and self.organism != 'All organisms':
                    if row[self.sbtab.columns_dict['!Organism']] != self.organism:
                        continue
                    if row[self.sbtab.columns_dict['!Organism']] != '':
                        continue
                    
                # exclude entries without a numeric value
                if row[mean_column] == '':
                    continue
                    
                # Michaelis constants need reaction AND species
                if row[self.sbtab.columns_dict['!QuantityType']] == 'Michaelis constant':
                    if row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']] not in self.reaction_list:
                        continue
                    elif row[self.sbtab.columns_dict['!Compound:SBML:species:id']] not in self.species_list:
                        continue

                # convert amount to concentration (for enzyme concentrations)
                if row[self.sbtab.columns_dict['!QuantityType']] == 'concentration of enzyme':
                    if row[self.sbtab.columns_dict['!Unit']] == 'molecules/cell':
                        continue
                        
                # exclude values that lie outside of the given boundaries;
                # this exclusion is triggered by flag ig_bounds
                if 'boundary_values' in self.parameter_dict and ig_bounds:
                    if self.parameter_dict['boundary_values'] == 'ignore':
                        if row[self.sbtab.columns_dict['!QuantityType']] != '':
                            if self.pmin[row[self.sbtab.columns_dict['!QuantityType']]] != None:
                                if float(row[mean_column]) < float(self.pmin[row[self.sbtab.columns_dict['!QuantityType']]]):
                                    if log_header == False:
                                        self.log += '\n### Warnings about ignored values that '\
                                                    'lie out of the boundaries ### \n'
                                        log_header = True
                                    self.log += 'The value %s for the %s of %s, %s lies below t'\
                                                'he requested minimum value of %s. It is ignor'\
                                                'ed for the balancing.\n'%(row[mean_column],
                                                                           row[self.sbtab.columns_dict['!QuantityType']],
                                                                           row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']],
                                                                           row[self.sbtab.columns_dict['!Compound:SBML:species:id']],
                                                                           self.pmin[row[self.sbtab.columns_dict['!QuantityType']]])
                                    continue
                        if row[self.sbtab.columns_dict['!QuantityType']] != '':
                            if self.pmax[row[self.sbtab.columns_dict['!QuantityType']]] != None:
                                if float(row[mean_column]) > float(self.pmax[row[self.sbtab.columns_dict['!QuantityType']]]):
                                    if log_header == False:
                                        self.log += '\n### Warnings about ignored values that '\
                                                    'lie out of the boundaries ### \n'
                                        log_header = True
                                    self.log += 'The value %s for the %s of %s, %s lies above t'\
                                                'he requested maximum value of %s. It is ignor'\
                                                'ed for the balancing.\n'%(row[mean_column],
                                                                           row[self.sbtab.columns_dict['!QuantityType']],
                                                                           row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']],
                                                                           row[self.sbtab.columns_dict['!Compound:SBML:species:id']],
                                                                           self.pmax[row[self.sbtab.columns_dict['!QuantityType']]])
                                    continue
                                
                # this part is only required when balancing with fixed
                # equilibrium constants and concentrations
                '''
                if row[self.sbtab.columns_dict['!QuantityType']] == 'concentration':
                    row[std_column] = float(0.00000001)
                    row[self.sbtab.columns_dict['!Min']] = row[mean_column]
                    row[self.sbtab.columns_dict['!Max']] = row[mean_column]
                elif row[self.sbtab.columns_dict['!QuantityType']] == 'equilibrium constant':
                    row[std_column] = float(0.00000001)
                    row[self.sbtab.columns_dict['!Min']] = row[mean_column]
                    row[self.sbtab.columns_dict['!Max']] = row[mean_column]
                elif row[self.sbtab.columns_dict['!QuantityType']] == 'reaction affinity':
                    row[std_column] = 0.01
                elif row[self.sbtab.columns_dict['!QuantityType']] == 'standard chemical potential':
                    row[std_column] = 1
                else:
                    row[self.sbtab.columns_dict['!Min']] = ''
                    row[self.sbtab.columns_dict['!Max']] = ''
                '''
                                  
                # reaction entities must not have assigned species, just like
                # species entities must not have assigned reaction
                if row[self.sbtab.columns_dict['!QuantityType']] in self.species_parameters:
                    row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']] = ''
                elif row[self.sbtab.columns_dict['!QuantityType']] in self.reaction_parameters:
                    row[self.sbtab.columns_dict['!Compound:SBML:species:id']] = ''

                # if the std is missing, use default std from prior file
                if row[std_column] == '':
                    row[std_column] = self.data_std[row[self.sbtab.columns_dict['!QuantityType']]]

                consistent_rows.append(row)
                
        return consistent_rows

    def new_row(self, name, quantity):
        '''
        generates one new row that is required for model structure
        '''
        row = ['']*len(self.new_header)
        row[0] = quantity
        if quantity in self.reaction_species_parameters:
            row[1] = name[1]
            row[2] = name[2]
        elif quantity in self.species_parameters:
            row[2] = name
        else: row[1] = name
        row[4] = self.quantity_type2unit[quantity]

        return row

    def existing_row(self, name, quantity):
        '''
        generates one row that is required and data are provided
        '''
        # if multiple entries are provided, they must be averaged
        value_dict = False
        amount = 1  # this is not even changed ever; is it required?
        if quantity in self.reaction_species_parameters:
            if self.available_parameters.count([quantity, name[1], name[2]]) > 1:
                value_dict = self.mean_row(name, quantity)
        elif self.available_parameters.count([quantity, name]) > 1:
            value_dict = self.mean_row(name, quantity)

        # after the averaging, the row is built
        used_rows = []
        new_row = [''] * len(self.new_header)

        for row in self.rows:
            for i in range(0, amount):
                further = False
                if quantity in self.reaction_species_parameters:
                    if row[self.sbtab.columns_dict['!QuantityType']] == quantity \
                       and row[self.sbtab.columns_dict['!Compound:SBML:species:id']] == name[2] \
                       and row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']] == name[1]:
                        new_row[1] = row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']]
                        new_row[2] = row[self.sbtab.columns_dict['!Compound:SBML:species:id']]
                        further    = True
                else:
                    if row[self.sbtab.columns_dict['!QuantityType']] == quantity \
                       and (row[self.sbtab.columns_dict['!Compound:SBML:species:id']] == name \
                            or row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']] == name) \
                            and not row in used_rows:
                        new_row[1] = row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']]
                        new_row[2] = row[self.sbtab.columns_dict['!Compound:SBML:species:id']]
                        further    = True
                if further:
                    new_row[0] = row[self.sbtab.columns_dict['!QuantityType']]
                    if value_dict:
                        new_row[3] = str(value_dict['Mode'])
                        new_row[5] = str(value_dict['Mean'])
                        new_row[6] = str(value_dict['Std'])
                    else:
                        new_row[5] = str(row[self.sbtab.columns_dict['!Mean']])
                        new_row[6] = str(row[self.sbtab.columns_dict['!Std']])
                        new_row[3] = str(round(self.normal_to_log([float(new_row[5])], [float(new_row[6])], [new_row[0]])[0][0], 4))

                    if quantity in self.thermodynamics: new_row[3] = new_row[5]
                    new_row[4] = self.quantity_type2unit[row[self.sbtab.columns_dict['!QuantityType']]]
                    # optional columns (columns 9 and more)

                    # MIN and MAX is currently out of order. Reinstall later.
                    if '!Min' in self.sbtab.columns_dict and '!Max' in self.sbtab.columns_dict:
                        new_row[j] = row[self.sbtab.columns_dict['!Min']]
                        new_row[j+1] = row[self.sbtab.columns_dict['!Max']]

                    # required?
                    if amount > 1 and i < amount:
                        multi_rows.append(new_row)
                        new_row = ['']*len(self.new_header)
                        used_rows.append(row)
                        if len(multi_rows) == amount:
                            return multi_rows

        return new_row

    def mean_row(self, name, quantity):
        '''
        if there are more than one value for one parameter, calculate the mean
        '''
        # collect available means and stds
        means = []
        stds = []
        for row in self.rows:
            if quantity in self.reaction_species_parameters:
                if row[self.sbtab.columns_dict['!QuantityType']] == quantity \
                   and row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']] == name[1] \
                   and row[self.sbtab.columns_dict['!Compound:SBML:species:id']] == name[2]:
                    means.append(float(row[self.sbtab.columns_dict['!Mean']]))
                    stds.append(float(row[self.sbtab.columns_dict['!Std']]))
            else:
                if row[self.sbtab.columns_dict['!QuantityType']] == quantity \
                   and (row[self.sbtab.columns_dict['!Compound:SBML:species:id']] == name \
                        or row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']] == name):
                    means.append(float(row[self.sbtab.columns_dict['!Mean']]))
                    stds.append(float(row[self.sbtab.columns_dict['!Std']]))
                    
        # build the mean
        if quantity in self.quantity_type2median_std:
            # geometric mean for all multiplicative quantities
            import scipy.stats
            denominator = 0
            for std in stds:
                if std != 0.0: denominator = denominator + (1/std**2)
                else: denominator = denominator + (1/numpy.log(2)**2)
            mean   = scipy.mean(means)
            std    = scipy.mean(stds)
            if not quantity in self.thermodynamics:
                median = numpy.exp(self.normal_to_log([mean], [std], False)[0])[0]
            else: median = mean
            value_dict = dict([('Mean', mean), ('Std', std), ('Mode', median)])
        else:
            # arithmetic mean for all additive/thermodynamic quantities
            denominator = 0
            for i, std in enumerate(stds):
                if std != 0.0: denominator = denominator + (1 / std ** 2)
                else: denominator = denominator + (1 / numpy.log(2) ** 2)

            std = numpy.sqrt(1/denominator)
            numerator   = 0
            denominator = 0
            for i, mean in enumerate(means):
                numerator   = numerator + mean / stds[i] ** 2
                denominator = denominator + 1 / stds[i] ** 2

            mean = numerator / denominator
            if not quantity in self.thermodynamics:
                median = numpy.exp(self.normal_to_log([mean], [std], False)[0])[0]
            else: median = mean
            value_dict = dict([('Mean', mean), ('Std', std), ('Mode', median)])
        
        return value_dict

    def fill_sbtab(self, sbtab, pseudos=None, priors=None):
        '''
        fills the values in the given SBtabfile
        '''
        self.pseudo_used = False
        self.make_default_table()
        sbtab_strings = [sbtab.header_row, '\t'.join(sbtab.columns)]
        
        
        if pseudos:
            # first fill parameter rows that have no value
            self.pseudo_used = True
            for i, row in enumerate(sbtab.value_rows):
                if row[sbtab.columns_dict['!QuantityType']] in self.pseudo_list:
                    try:
                        row[sbtab.columns_dict['!Mode']] = str(pseudos[row[0]][0])
                        row[sbtab.columns_dict['!UnconstrainedGeometricMean']] = str(pseudos[row[0]][0])
                        row[sbtab.columns_dict['!UnconstrainedGeometricStd']] = str(pseudos[row[0]][1])
                        sbtab_strings.append('\t'.join(row))
                    except: pass
                    
            # then construct required variables
            means = []
            stds = []
            for i, quantity in enumerate(self.pseudo_list):
                means.append(pseudos[quantity][0])
                stds.append(pseudos[quantity][1])
                
            (self.log_means, self.log_stds) = self.med10_std_to_log(means, stds, self.pseudo_list)
             
            sbtab_pseudo = SBtab.SBtabTable('\n'.join(sbtab_strings), 'sbtab_pseudo.csv')
            return sbtab_pseudo
                           
        return sbtab

    def make_default_table(self):
        '''
        generate a table of values for every parameter for every reaction/species
        '''
        medians = []
        stds = []
        quantities = []
        counter = 0
        
        for parameter_type in self.thermodynamics:
            medians.append(self.quantity_type2mean_std[parameter_type][0])
            quantities.append(parameter_type)
            stds.append(self.quantity_type2mean_std[parameter_type][1])
            counter += 1
            
        for parameter_type in self.quantity_type2median_std.keys():
            medians.append(self.quantity_type2median_std[parameter_type][0])
            stds.append(self.quantity_type2median_std[parameter_type][1])
            quantities.append(parameter_type)
            counter += 1
            
        (self.log_means, self.log_stds) = self.med10_std_to_log(medians, stds, quantities)
        
        self.prior_values = {}
        for i, quantity in enumerate(quantities):
            self.prior_values[quantity] = [(self.log_means[i], self.log_stds[i])]
            
    def normal_to_log(self, means, stds, types):
        '''
        generates log values for normal values
        '''
        log_means = []
        log_stds  = []

        for i, mean in enumerate(means):
            if types and types[i] in self.thermodynamics:
                log_means.append(float(mean))
                log_stds.append(float(stds[i]))
                
            else:
                term = numpy.log(1+(numpy.square(float(stds[i]))/numpy.square(float(mean))))
                log_means.append(numpy.log(float(mean))-0.5*term)
                log_stds.append(numpy.sqrt(numpy.log(1+(numpy.square(float(stds[i]))/numpy.square(float(mean))))))
               
        if 'nan' in log_means:
            raise ParameterBalancingError('The logarithm of one of your given mean values is invalid.')

        return log_means, log_stds

    def log_to_normal(self, log_means, log_stds, types=None):
        '''
        generates a list of the normal values from a list of the log values
        '''
        means = []
        stds  = []
        for i, log_mean in enumerate(log_means):
            if types and types[i] in self.thermodynamics:
                means.append(log_mean)
                stds.append(log_stds[i])
            else:
                means.append(numpy.exp(log_mean+0.5*numpy.square(log_stds[i])))
                stds.append(numpy.sqrt((numpy.exp(numpy.square(float(log_stds[i])))-1) * (numpy.exp(2*log_mean+numpy.square(float(log_stds[i]))))))
        
        return means, stds

    def med10_std_to_log(self, medians, stdlogs, types):
        '''
        this is the users default choice for missing values: inserted are the median and
        the stdlog10, the output are the corresponding mean value and the standard dev.
        '''
        log_means = []
        log_stds  = []

        for i, median in enumerate(medians):
            if types and types[i] in self.thermodynamics:
                log_means.append(median)
                log_stds.append(stdlogs[i])
            else:
                log_means.append(numpy.log(median))
                log_stds.append(numpy.log(stdlogs[i]))

        return log_means, log_stds

    def add_config_to_log(self):
        '''
        adds config information for the balancing
        '''
        self.log += '\n### You have used an options file for the parameter balancing. ### \nThe options are as follows:\n'
        self.log += '!Option\t!Value\n'
        for entry in self.parameter_dict.keys():
            self.log += '%s\t%s\n'%(entry, self.parameter_dict[entry])

###############################################################################

    def print_warning(self, type, flag):
        if self.warned:
            print('There was an error in the numerics. This may be due to broadly chosen probability distributions. Try the usage of pseudo values for a fix.')
            self.warned = False

    def make_balancing(self, sbtab, sbtab_old, pmin, pmax, parameter_dict):
        '''
        generates the values for the parameter balancing
        '''
        self.sbtab = sbtab_old
        self.sbtab_new = sbtab

        # initialise log file
        self.parameter_dict = parameter_dict
        if 'config' in self.parameter_dict.keys():
            self.add_config_to_log()
            
        # initialise needed variables
        self.pmin             = pmin
        self.pmax             = pmax
        self.new_header       = header_names

        new_rows = self.sbtab_new.value_rows
        self.new_rows = self.sort_list(new_rows)

        # get matrix information from hardcoded sheet
        self.sheet            = self.get_sheet()

        self.temperature      = float(self.parameter_dict['temperature'])
        self.pH               = float(self.parameter_dict['ph'])

        # build needed vectors and matrices
        self.desired_parameters = self.build_desired_parameters()
        self.x_vector           = self.collect_available_values()
        self.theta_vector       = self.build_theta_vector()
        (self.C_prior, self.C_x) = self.build_covariance_matrices()
        
        # build the dependence matrix D and the data specific D_x
        self.Q   = self.build_dependence_matrix()
        self.Q_star = self.build_specific_dependence_matrix()

        # KEY PART: calculating the posteriori covariance matrix C_post and the posteriori mean vector mean_post
        self.calculate_posteriori()

        # make normal values again
        (self.mean_post, self.stds_post) = self.log_to_normal(self.x_post, self.stds_log_post, self.quantities)

        ################################################################
        # generating minimization problem
        self.optimized = False
        if self.bounds_inc.count((None, None)) != len(self.bounds_inc):  # and False:    # and False ---> mustn't be loaded in web interface
            # generating the medians to bound them                                       # because the files have to be generated manually!!
            print('in optimiser')                                                        # only for super-cool users...
            
            medians = []
            medstds = []

            (self.means_inc, self.stds_inc) = self.log_to_normal(self.means_post_inc, self.stds_log_inc, self.quantities_inc)

            for i, value in enumerate(self.means_inc):
                if not self.quantities_inc[i] in self.thermodynamics:
                    (log_mean, log_std) = self.normal_to_log([self.means_inc[i]], [self.stds_inc[i]], self.quantities_inc[i])
                    medians.append(numpy.exp(log_mean[0]))
                else:
                    medians.append(self.means_inc[i])
                medstds.append(max(self.stds_inc[i], 10))

            # optimisation function
            def log_mean_post_func(q):
                return numpy.dot(((q-medians).transpose()), (numpy.dot(numpy.linalg.inv(self.C_post), (q-medians))))

            # setting proper boundaries for parameters that have no boundaries set in SBtab
            new_boundaries = []
            is_logarithmic = []
            for i, bound in enumerate(self.bounds_inc):
                if bound == ('', '') or bound == (None, None):
                    # setting boundaries for thermodynamic parameters (kJ/mol)
                    # THIS NEEDS TO BE REWRITTEN DESPERATELY;
                    # NO USAGE OF NAME2INDEX OR THERMODYNAMIC_INDICES;
                    if name2index[self.quantities_inc[i]] in thermodynamic_indices:
                        is_logarithmic.append(False)
                        if (max(medians[i]-medstds[i]*2, -3000))<(min(3000, medians[i]+medstds[i]*2)):
                            new_boundaries.append((max(medians[i]-medstds[i]*2, -3000), min(3000, medians[i]+medstds[i]*2)))
                        else:
                            new_boundaries.append((medians[i]+medstds[i]*2, medians[i]-medstds[i]*2))
                    # setting boundaries for all other parameters
                    else:
                        is_logarithmic.append(True)
                        if (medians[i]-medstds[i]*4)<(medians[i]+medstds[i]*4):
                            new_boundaries.append((max(medians[i]-medstds[i]*4, 0.00001), medians[i]+medstds[i]*4))
                        else:
                            new_boundaries.append((medians[i]+medstds[i]*4, medians[i]-medstds[i]*4))
                else:
                    is_logarithmic.append(True)
                    new_boundaries.append((float(bound[0]), float(bound[1])))

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

            f = open('medians.txt', 'w')
            for i, element in enumerate(medians):
                if not i == len(medians)-1:
                    f.write(str(element)+',')
                else:
                    f.write(str(element))
            f.close()

            g = open('cpost.txt', 'w')
            for line in self.C_post:
                for i, element in enumerate(line):
                    if not i == len(line)-1:
                        g.write(str(element)+',')
                    else:
                        g.write(str(element))
                g.write('\n')
            g.close()

            # generating optimization and updating responsible mean vector
            new_medians = misc.fmin_gen(log_mean_post_func, numpy.array(medians), population_size=20, survivors=5, generations=500, bounds=proper_boundaries, use_pp=False, variable_is_logarithmic=is_logarithmic, disp=1)
            for i, value in enumerate(new_medians):
                if value+0.001 < proper_boundaries[i][0] or value-0.001 > proper_boundaries[i][1]:
                    self.quantities_inc[i]
                    print('NEW_MODE value out of bound: ', value, ' [bounds: ', proper_boundaries[i], ']')

            (new_medians_log, new_stds_log) = self.normal_to_log(new_medians, self.stds_inc, self.quantities_inc)
            self.C_xpost = numpy.dot((numpy.dot(self.Q, self.C_post)), self.Q.transpose())
            self.x_post = numpy.dot(self.Q, new_medians_log)
            self.stds_log_post = self.extract_cpost()

            (self.mean_post_opt, self.stds_post_opt) = self.log_to_normal(self.x_post, self.stds_log_post, self.quantities)
            
            self.optimized = True
      
        #################################################################
        # make value-dictionaries to realize the insertion of the computed values into the SBtab-file in GUI
        balanced_sbtab = self.build_new_sbtab()

        sbtab_string = []
        for entry in balanced_sbtab:
            sbtab_string.append('\t'.join(entry))
            
        sbtab_new = SBtab.SBtabTable('\n'.join(sbtab_string), 'sbtab_new.csv')
        C_string = self.make_cpost_string()
        shannons = self.get_shannons()

        return sbtab_new, self.mean_post, self.q_post, C_string, self.C_post, self.Q, shannons, self.log

    def get_sheet(self):
        '''
        get the sheet that tells us, how to build up which matrix and further stuff
        '''
        sheet = {"equilibrium constant": ((-1/2.4790, "A"), 0, 0, 0, 0, 0, 0),
                 "substrate catalytic rate constant": ((-0.5/2.4790, "A"), 1, (-0.5, "Z"), 0, 0, 0, 0),
                 "product catalytic rate constant": ((0.5/2.4790, "A"), 1, (0.5, "Z"), 0, 0, 0, 0),
                 "forward maximal velocity": ((-0.5/2.4790, "A"), 1, (-0.5, "Z"), 0, 0, 1, 0),
                 "reverse maximal velocity": ((0.5/2.4790, "A"), 1, (0.5, "Z"), 0, 0, 1, 0),
                 "chemical potential": (1, 0, 0, 0, 0, 0, (2.790, '1')),
                 "reaction affinity": ((-1, "A"), 0, 0, 0, 0, 0, (-2.4790, "AB"))}

        return sheet

    def build_desired_parameters(self):
        '''
        builds up a dictionary of identifiers for those parameters that are chosen to
        be used for the balancing.
        they are of the kind: [(equilibrium constant, reaction)]=keq_1
        '''
        desired_parameters       = {}
        for parameter_type in self.parameter_dict.keys():
            if self.parameter_dict[parameter_type] == True:
                if parameter_type in self.species_parameters:
                    for species in self.species_list:
                        desired_parameters[parameter_type, species] = self.quantity2identifier[parameter_type]+'_'+str(self.species2number[species])
                elif parameter_type in self.reaction_parameters:
                    for reaction in self.reaction_list:
                        desired_parameters[parameter_type, reaction] = self.quantity2identifier[parameter_type]+'_'+str(self.reaction2number[reaction])
                elif parameter_type in self.reaction_species_parameters:
                    for reaction_species in self.model_specific:
                        if parameter_type == reaction_species[0]:
                            desired_parameters[parameter_type, reaction_species[1], reaction_species[2]] = self.quantity2identifier[parameter_type]+'_'+str(self.reaction2number[reaction_species[1]])+'_'+str(self.species2number[reaction_species[2]])

        return desired_parameters

    def collect_available_values(self):
        '''
        generates the x-vector (values that are provided by the SBtab).
        we have to get the values from the initial SBtab first, afterwards
        the values that have been added by the use of priors and pseudos.
        so we are basically searching two SBtab files and storing their content.
        '''
        self.quantities_x = []
        self.parameter2bounds = {}

        # first: get the column indices of the initial SBtab file, so we can search it
        old_rows = self.tidy_up_sbtab(True)

        value_tuples_old = []
        value_tuples_new = []
        row_identifiers = []

        try: mean_column = self.sbtab.columns_dict['!Mean']
        except: mean_column = self.sbtab.columns_dict['!UnconstrainedGeometricMean']
        try: std_column = self.sbtab.columns_dict['!Std']
        except: std_column = self.sbtab.columns_dict['!UnconstrainedGeometricStd']

        for row in old_rows:
            quantity = row[self.sbtab.columns_dict['!QuantityType']]
            if quantity in self.parameter_dict.keys():
                if len(row) == len(self.sbtab.columns):
                    single_tuple = []
                    single_tuple.append(row[self.sbtab.columns_dict['!QuantityType']])
                    single_tuple.append(row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']])
                    single_tuple.append(row[self.sbtab.columns_dict['!Compound:SBML:species:id']])
                    single_tuple.append(row[mean_column])
                    if not row[std_column] == '': single_tuple.append(row[std_column])
                    elif row[self.sbtab.columns_dict['!QuantityType']] in self.thermodynamics: single_tuple.append('35.0')
                    else: single_tuple.append(str(float(row[mean_column])*0.5))
                    single_tuple.append(self.make_identifier(row))
                    self.quantities_x.append(row[self.sbtab.columns_dict['!QuantityType']])
                   
                    # check, whether this parameter is really part of the model structure and store boundaries
                    if not row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']] == '' and not row[self.sbtab.columns_dict['!Compound:SBML:species:id']] == '':
                        row_identifier = (row[self.sbtab.columns_dict['!QuantityType']], row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']], row[self.sbtab.columns_dict['!Compound:SBML:species:id']])
                        try:
                            if self.sbtab.columns_dict['!Min'] and self.sbtab.columns_dict['!Max']:
                                self.parameter2bounds[row[self.sbtab.columns_dict['!QuantityType']], (row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']], row[self.sbtab.columns_dict['!Compound:SBML:species:id']])] = (row[self.sbtab.columns_dict['!Min']], row[self.sbtab.columns_dict['!Max']])
                            else: self.parameter2bounds[row[self.sbtab.columns_dict['!QuantityType']], (row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']], row[self.sbtab.columns_dict['!Compound:SBML:species:id']])] = (None, None)
                        except: pass
                    elif not row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']] == '':
                        row_identifier = (row[self.sbtab.columns_dict['!QuantityType']], row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']])
                        try:
                            if self.sbtab.columns_dict['!Min'] and self.sbtab.columns_dict['!Max']:
                                self.parameter2bounds[row[self.sbtab.columns_dict['!QuantityType']], row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']]] = (row[self.sbtab.columns_dict['!Min']], row[self.sbtab.columns_dict['!Max']])
                            else: self.parameter2bounds[row[self.sbtab.columns_dict['!QuantityType']], (row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']])] = (None, None)
                        except: pass
                    else:
                        row_identifier = (row[self.sbtab.columns_dict['!QuantityType']], row[self.sbtab.columns_dict['!Compound:SBML:species:id']])
                        try:
                            if self.sbtab.columns_dict['!Min'] and self.sbtab.columns_dict['!Max']:
                                self.parameter2bounds[row[self.sbtab.columns_dict['!QuantityType']], row[self.sbtab.columns_dict['!Compound:SBML:species:id']]] = (row[self.sbtab.columns_dict['!Min']], row[self.sbtab.columns_dict['!Max']])
                            else: self.parameter2bounds[row[self.sbtab.columns_dict['!QuantityType']], (row[self.sbtab.columns_dict['!Compound:SBML:species:id']])] = (None, None)
                        except: pass
                    row_identifiers.append(row_identifier)
                    if row_identifier in self.desired_parameters.keys(): value_tuples_old.append(single_tuple)
            else:
                raise ParameterBalancingError('There is a quantity type in the SBtab that cannot be interpreted: %s'%(row[self.sbtab.columns_dict['!QuantityType']]))

        # second: get column indices of the new SBtab file and search it
        for row in self.new_rows:
            if row[self.sbtab_new.columns_dict['!UnconstrainedGeometricMean']] == '': continue
            if len(row) == len(self.sbtab_new.columns):
                single_tuple = []
                single_tuple.append(row[self.sbtab_new.columns_dict['!QuantityType']])
                single_tuple.append(row[self.sbtab_new.columns_dict['!Reaction:SBML:reaction:id']])
                single_tuple.append(row[self.sbtab_new.columns_dict['!Compound:SBML:species:id']])
                if not row[self.sbtab_new.columns_dict['!Reaction:SBML:reaction:id']] == '' and not row[self.sbtab_new.columns_dict['!Compound:SBML:species:id']] == '':
                    row_identifier = (row[self.sbtab_new.columns_dict['!QuantityType']],
                                      (row[self.sbtab_new.columns_dict['!Reaction:SBML:reaction:id']],
                                       row[self.sbtab_new.columns_dict['!Compound:SBML:species:id']]))
                elif not row[self.sbtab_new.columns_dict['!Reaction:SBML:reaction:id']] == '':
                    row_identifier = (row[self.sbtab_new.columns_dict['!QuantityType']], row[self.sbtab_new.columns_dict['!Reaction:SBML:reaction:id']])
                else: row_identifier = (row[self.sbtab_new.columns_dict['!QuantityType']], row[self.sbtab_new.columns_dict['!Compound:SBML:species:id']])

                if not row_identifier in row_identifiers:
                    try:
                        if self.sbtab_new.columns_dict['!Min'] and self.sbtab_new.columns_dict['!Max']:
                            self.parameter2bounds[row_identifier] = (None, None)
                    except: pass
                if row[self.sbtab_new.columns_dict['!QuantityType']] in self.prior_list and not self.pseudo_used: continue
                elif row[self.sbtab_new.columns_dict['!QuantityType']] in self.pseudo_list and not self.pseudo_used: continue
                
                if not row_identifier in row_identifiers:
                    single_tuple.append(row[self.sbtab_new.columns_dict['!UnconstrainedGeometricMean']])
                    single_tuple.append(row[self.sbtab_new.columns_dict['!UnconstrainedGeometricStd']])
                    single_tuple.append(self.make_identifier(row))
                    try:
                        if self.sbtab_new.columns_dict['!Min'] and row[self.sbtab_new.columns_dict['!Min']] != '':
                            single_tuple.append((float(row[self.sbtab_new.columns_dict['!Min']]), float(row[self.sbtab_new.columns_dict['!Max']])))
                        elif self.sbtab_new.columns_dict['!Min']:
                            single_tuple.append((None, None))
                        if self.sbtab_new.columns_dict['!Min'] and self.sbtab_new.columns_dict['!Max']:
                            self.parameter2bounds[row_identifier] = (None, None)
                    except: pass
                    value_tuples_new.append(single_tuple)

        # generate logarithms
        means = []
        stds = []
        types = []
        vt = []
        self.x_star = []
       
        for single_tuple in value_tuples_old:
            if single_tuple[3] == '': continue
            means.append(float(single_tuple[3]))
            stds.append(float(single_tuple[4]))
            types.append(single_tuple[0])
            vt.append(single_tuple)
            (self.x_star, self.log_stds_x) = self.normal_to_log(means, stds, types)

        self.new_rows = self.sbtab.value_rows + self.new_rows
        
        return vt

    def sort_list(self, qlist):

        quantities = ['standard chemical potential',
                      'catalytic rate constant geometric mean',
                      'Michaelis constant', 'activation constant',
                      'inhibitory constant', 'concentration of enzyme',
                      'concentration', 'equilibrium constant',
                      'substrate catalytic rate constant',
                      'product catalytic rate constant',
                      'forward maximal velocity',
                      'reverse maximal velocity', 'chemical potential',
                      'reaction affinity']

        new_list = []
        for quantity in quantities:
            for entry in qlist:
                if type(entry) == list:
                    if entry[0] == quantity:
                        new_list.append(entry)
                else:
                    if entry == quantity:
                        new_list.append(quantity)
        return new_list

    def make_identifier(self, row):
        '''
        generates a specific identifier that works for only this specific row
        '''
        # qidentifier: kG, kV, kM, etc.
        qidentifier = self.quantity2identifier[row[self.sbtab.columns_dict['!QuantityType']]]
        # qnumber: specific species and/or reaction number in the order of their appearance in SBML model
        if row[self.sbtab.columns_dict['!QuantityType']] in self.species_parameters: qnumber = str(self.species2number[row[self.sbtab.columns_dict['!Compound:SBML:species:id']]])
        elif row[self.sbtab.columns_dict['!QuantityType']] in self.reaction_parameters: qnumber = str(self.reaction2number[row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']]])
        else: qnumber = str(self.species2number[row[self.sbtab.columns_dict['!Compound:SBML:species:id']]])+'_'+str(self.reaction2number[row[self.sbtab.columns_dict['!Reaction:SBML:reaction:id']]])
        identifier = qidentifier+'_'+qnumber

        return identifier

    def get_shannons(self):
        '''
        calculate the shannon entropies of the prior and the posterior covariance matrices
        '''
        shannons = []
        # first, shannon entropy of prior matrix self.C_prior
        shannon_prior = 0.5 * numpy.log((2*numpy.pi*numpy.exp(2)) * scipy.linalg.det(self.C_prior))
        shannons.append(shannon_prior)

        # second, get shannon entropy of posterior matrix
        shannon_posterior = 0.5 * numpy.log((2*numpy.pi*numpy.exp(2)) * scipy.linalg.det(self.C_post))
        shannons.append(shannon_posterior)

        return shannons

    def make_cpost_string(self):
        '''
        in order to export the posterior covariance matrix to the user, we make a tsv-string out of it (nicer to read)
        '''
        new_C = '\n'.join(('\t'.join(str(e) for e in row)) for row in self.C_xpost)

        return new_C

    def sample_posterior_dist(self, posterior, posterior_inc, c_matrix_inc, r_matrix, new_SBtab, header, number):
        '''
        give the posterior to this function and receive *number* times a sample from the posterior distribution
        '''
        # first, get the matrix root of the posterior covariance matrix
        C_root = scipy.linalg.matfuncs.sqrtm(c_matrix_inc)
        list_of_SBtab_strings = []

        # do for *number* times:
        for i in range(number):
            # generate random variables of posteriors length
            ksi = numpy.random.normal(0, 1, len(posterior_inc))

            # third: get the new posterior for the new SBtab
            new_posterior = numpy.array(posterior)+(numpy.dot(numpy.dot(r_matrix, C_root), ksi))

            # last, but not least: get the new SBtab with the sampled values
            new_SBtab = self.build_new_sbtab(new_sbtab=new_SBtab, posterior_sample=new_posterior, header=header)

            list_of_SBtab_strings.append(new_SBtab)

        return list_of_SBtab_strings 

    def build_dependence_matrix(self):
        '''
        builds the dependence matrix D from all the submatrices needed for the balancing
        '''
        # the dependence matrix consists of two major parts: the unit matrix on top
        # and the rows below the unit matrix
        D_matrix = []
        self.parameter2row = {}
        self.quantities = []
        
        # first, we build up the unit matrix
        unit_matrix = self.build_unit_matrix()
        for row in unit_matrix:
            D_matrix.append(row)

        # second, we check which bottom rows we have to build up
        for pseudo_quantity in self.pseudo_list:
            if self.parameter_dict[pseudo_quantity]:
                rows = self.build_bottom_row(pseudo_quantity)
                for row in rows:
                    D_matrix.append(row)

        matrix = numpy.array(D_matrix)
        return matrix

    def build_unit_matrix(self):
        '''
        builds up the unit matrix as a first part of the dependence matrix D
        uses only the prior parameters that are chosen by the user
        '''
        unit_rows       = []
        self.bounds     = []
        self.id_order   = {}

        for i, x in enumerate(self.theta_basic):
            row    = [0.0]*len(self.theta_basic)
            row[i] = 1.0
            unit_rows.append(row)
            self.parameter2row[(x[0], x[2])] = row
            self.quantities.append(x[0])
            self.id_order[(x[0], x[2])] = i
            if '!Min' in self.sbtab_new.columns_dict:
                self.bounds.append(self.parameter2bounds[(x[0], x[2])])
                
        self.matrix_row_counter = len(self.id_order)
        return unit_rows

    def create_row_specifics(self, row_specifics_str):
        '''
        create a computable list of row specifics from the given string
        inherited from the prior file
        '''
        variable2value = {'RT': '2.4790', 'Nt': 'A'}        
    
    def build_bottom_row(self, pseudo_quantity):
        '''
        builds one of the bottom rows of the dependence matrix D individually.
        of course, this function is not short. but the build up of the specific
        rows of D is extremely complex and shall be able to be generic.
        '''
        #######################
        # building matrix info dynamically
        # let's try this later
        row_specifics_str = self.matrix_info[pseudo_quantity]
        row_specifics = self.create_row_specifics(row_specifics_str)
        #######################
        
        use_list = self.quantity2list[pseudo_quantity]       
        sheet = self.sheet[pseudo_quantity]
        self.remember_links = {}
        rows = []

        row_index = 0
        for i, element in enumerate(use_list):
            row          = [0.0]*len(self.theta_basic)
            column_index = 0
            if '!Min' in self.sbtab.columns_dict and '!Max' in self.sbtab.columns_dict:
                self.bounds.append(self.parameter2bounds[pseudo_quantity, element])
            for j, matrix_type in enumerate(sheet):
                if self.parameter_dict[self.prior_list[j]]:
                    index_list = self.quantity2list[self.prior_list[j]]

                    # build zero matrix
                    if matrix_type == 0:
                        column_index += (len(index_list))

                    # build unit matrix
                    if matrix_type == 1:
                        for k in range(len(index_list)):
                            row[column_index+row_index] = 1.0
                            column_index += 1
                            break
                        column_index += len(index_list)-1

                    # build specific matrices
                    try:
                        if len(matrix_type)>1:
                            factor = matrix_type[0]  # this factor represents e.g. R*T or 1/R*T
                            matrix = matrix_type[1]
                            try:
                                reactants       = self.reactions_reactants[element][0]
                                stoichiometry_r = self.reactions_reactants[element][1]
                                products        = self.reactions_products[element][0]
                                stoichiometry_p = self.reactions_products[element][1]
                            except: pass

                            # build N, the stoichiometric matrix
                            if matrix == 'A':
                                for species in self.species_list:
                                    for r, reactant in enumerate(reactants):
                                        if species == reactant:
                                            row[column_index] = stoichiometry_r[r]*factor
                                            column_index += 1
                                    for p, product in enumerate(products):
                                        if species == product:
                                            row[column_index] = stoichiometry_p[p]*factor
                                            column_index += 1

                            # for the end of reaction affinities
                            elif matrix == 'AB':
                                for species in self.species_list:
                                    for r, reactant in enumerate(reactants):
                                        if species == reactant:
                                            row[column_index] = stoichiometry_r[r]*factor
                                            column_index += 1
                                    for p, product in enumerate(products):
                                        if species == product:
                                            row[column_index] = stoichiometry_p[p]*factor
                                            column_index += 1

                            # build Z, the values for the Michaelis constant coefficients
                            elif matrix == 'Z':
                                if matrix_type[0] < 0: factor = -1.0
                                else: factor = 1.0
                                for michaelis_tuple in self.model_michaelis:
                                    if michaelis_tuple[2] in reactants:
                                        row[column_index] = -0.5*factor
                                        column_index += 1
                                    if michaelis_tuple[2] in products:
                                        row[column_index] = 0.5*factor
                                        column_index += 1
                                    if michaelis_tuple[2] not in reactants and michaelis_tuple[2] not in products:
                                        column_index += 1

                            # build 1, a simple alternative to N (for chemical potentials)
                            elif matrix == '1':
                                for k in range(len(index_list)):
                                    row[column_index+row_index] = factor
                                    column_index += 1
                                    break
                                column_index += len(index_list)-1
                    except: pass

            row_index += 1
            rows.append(row)
            self.quantities.append(pseudo_quantity)
            self.parameter2row[(pseudo_quantity, use_list[i])] = row
            self.id_order[(pseudo_quantity, element)] = self.matrix_row_counter
            self.matrix_row_counter += 1
            
        return rows

    def build_specific_dependence_matrix(self):
        '''
        the matrix D_x is the dependence matrix that holds all rows from D which can
        be also found in the x_vector. Thus, the dependence matrix for the available
        values.
        '''
        rows = []
        
        for i, single_tuple in enumerate(self.x_vector):
            if single_tuple[0] in self.species_parameters:
                rows.append(self.parameter2row[(single_tuple[0], single_tuple[2])])
            elif single_tuple[0] in self.reaction_parameters:
                rows.append(self.parameter2row[(single_tuple[0], single_tuple[1])])
            elif single_tuple[0] in self.reaction_species_parameters:
                rows.append(self.parameter2row[(single_tuple[0], (single_tuple[1], single_tuple[2]))])
            else: print('row identifier not found: ', single_tuple[0], ' ', single_tuple[1], ' ', single_tuple[2])

        if rows == []: matrix = 0
        else: matrix = numpy.array(rows)

        return matrix

    def build_theta_vector(self):
        '''
        generates the theta_vector (default prior means for every parameter in model)
        '''
        self.parameter2row = {}
        theta = []
        self.theta_basic = []
        self.q_prior = []
        self.log_stds_prior  = []
        types = []
        used_identifiers = []
        self.quantities_inc   = []
        self.bounds_inc  = []

        for quantity in self.prior_list:
            if self.parameter_dict[quantity]:
                if quantity in self.species_parameters:
                    for species in self.species_list:
                        if not (quantity, species) in used_identifiers:
                            theta.append((quantity, self.prior_values[quantity][0][0], species))
                            self.theta_basic.append((quantity, self.prior_values[quantity][0][0], species))
                            self.q_prior.append(self.prior_values[quantity][0][0])
                            self.log_stds_prior.append(self.prior_values[quantity][0][1])
                            types.append(quantity)
                            self.quantities_inc.append(quantity)
                            if (quantity, species) in self.parameter2bounds.keys():
                                if self.min_column:
                                    self.bounds_inc.append(self.parameter2bounds[(quantity, species)])
                            used_identifiers.append((quantity, species))
                elif quantity in self.reaction_parameters:
                    for reaction in self.reaction_list:
                        if not (quantity, reaction) in used_identifiers:
                            theta.append((quantity, self.prior_values[quantity][0][0], reaction))
                            self.theta_basic.append((quantity, self.prior_values[quantity][0][0], reaction))
                            self.q_prior.append(self.prior_values[quantity][0][0])
                            self.log_stds_prior.append(self.prior_values[quantity][0][1])
                            types.append(quantity)
                            self.quantities_inc.append(quantity)
                            if (quantity, reaction) in self.parameter2bounds.keys():
                                if self.min_column:
                                    self.bounds_inc.append(self.parameter2bounds[(quantity, reaction)])
                            used_identifiers.append((quantity, reaction))
                elif quantity in self.reaction_species_parameters:
                    for reaction_species in self.model_dict[quantity]:
                        if not (quantity, (reaction_species[1], reaction_species[2])) in used_identifiers:                        
                            theta.append((quantity, self.prior_values[quantity][0][0], (reaction_species[1], reaction_species[2])))
                            self.theta_basic.append((quantity, self.prior_values[quantity][0][0], (reaction_species[1], reaction_species[2])))
                            self.q_prior.append(self.prior_values[quantity][0][0])
                            self.log_stds_prior.append(self.prior_values[quantity][0][1])
                            types.append(quantity)
                            self.quantities_inc.append(quantity)
                            if (quantity, (reaction_species[1], reaction_species[2])) in self.parameter2bounds.keys():
                                if self.min_column:
                                    self.bounds_inc.append(self.parameter2bounds[(quantity, (reaction_species[1], reaction_species[2]))])
                            used_identifiers.append((quantity, (reaction_species[1], reaction_species[2])))

        if self.pseudo_used:
            for quantity in self.pseudo_list:
                if self.parameter_dict[quantity]:
                    if quantity in self.species_parameters:
                        for species in self.species_list:
                            if not (quantity, species) in used_identifiers:
                                theta.append((quantity, self.prior_values[quantity][0][0], species))
                                self.quantities_inc.append(quantity)
                                if (quantity, species) in self.parameter2bounds.keys():
                                    if self.min_column:
                                        self.bounds_inc.append(self.parameter2bounds[(quantity, species)])
                                used_identifiers.append((quantity, species))
                                self.q_prior.append(self.prior_values[quantity][0][0])
                    elif quantity in self.reaction_parameters:
                        for reaction in self.reaction_list:
                            if not (quantity, reaction) in used_identifiers:
                                theta.append((quantity, self.prior_values[quantity][0][0], reaction))
                                self.quantities_inc.append(quantity)
                                if (quantity, reaction) in self.parameter2bounds.keys():
                                    if self.min_column:
                                        self.bounds_inc.append(self.parameter2bounds[(quantity, reaction)])
                                used_identifiers.append((quantity, reaction))
                                self.q_prior.append(self.prior_values[quantity][0][0])
                    elif quantity in self.reaction_species_parameters:
                        for reaction_species in self.model_dict[quantity]:
                            if not (quantity, (reaction_species[1], reaction_species[2])) in used_identifiers:
                                theta.append((quantity, self.prior_values[quantity][0][0], (reaction_species[1], reaction_species[2])))
                                self.quantities_inc.append(quantity)
                                if (quantity, (reaction_species[1], reaction_species[2])) in self.parameter2bounds.keys():
                                    if self.min_column:
                                        self.bounds_inc.append(self.parameter2bounds[(quantity, (reaction_species[1], reaction_species[2]))])
                                used_identifiers.append((quantity, (reaction_species[1], reaction_species[2])))
                                self.q_prior.append(self.prior_values[quantity][0][0])
        
        return theta

    def build_covariance_matrices(self):
        '''
        generate covariance matrix for measured values x (stds) from SBtab
        '''
        # first, generate prior covariance matrix C_prior
        
        C_prior_rows = []
        
        for i, theta in enumerate(self.theta_vector):
            row = [0.0]*len(self.theta_vector)
            row[i] = numpy.square(float(self.get_default_std(self.theta_vector[i][0])))
            C_prior_rows.append(row)

        C_prior = numpy.array(C_prior_rows)
        # second, generate covariance matrix according to the input values in the x-vector
        C_x_rows = []
        for i, x_entry in enumerate(self.x_vector):
            row = [0.0]*len(self.x_vector)
            row[i] = numpy.square(self.log_stds_x[i])
            C_x_rows.append(row)

        if C_x_rows == []:
            C_x = 0
        else:
            C_x = numpy.array(C_x_rows)

        return C_prior, C_x


    def get_default_std(self, quantity):
        '''
        returns the default standard deviation
        '''
        return self.prior_values[quantity][0][1]
    
    def calculate_posteriori(self):
        '''
        calculates the posteriori values
        '''
        # if no data is given, these variables are zero
        if self.x_star == []: self.x_star = 0

        try: self.C_x_inv = numpy.linalg.inv(self.C_x)
        except: self.C_x_inv = 0

        try: Q_star_trans = self.Q_star.transpose()
        except: Q_star_trans = 0
        
        # matrix inverse
        try: self.C_prior_inv = numpy.linalg.inv(self.C_prior)
        except:
            print("C_prior is not invertible\n")
            sys.exit()

        # for i, row in enumerate(self.C_prior):
        #    print(self.quantities_inc[i], ',', list(row))

        # for i, elem in enumerate(self.q_prior):
        #    print(self.quantities_inc[i], ',', elem)

        # for i, row in enumerate(self.Q):
        #    print(self.quantities[i], ',', list(row))

        # print(self.C_x)
        # for i, row in enumerate(self.quantities_x):
        #    print(row, ',', list(self.C_x[i]))

        # print(self.x_star)
        # for i, row in enumerate(self.x_star):
        #    print(self.quantities_x[i],  ',',  row)

        # print(self.Q_star)
        # for i, row in enumerate(self.Q_star):
        #    print(self.quantities_x[i],  ',', list(row))  
        
        # posterior covariance matrices
        if self.pseudo_used:
            self.C_post = numpy.linalg.inv(numpy.dot(numpy.dot(self.Q.transpose(), self.C_prior_inv), self.Q) + numpy.dot(numpy.dot(Q_star_trans, self.C_x_inv), self.Q_star))
        else:
            self.C_post = numpy.linalg.inv(self.C_prior_inv+numpy.dot(numpy.dot(Q_star_trans, self.C_x_inv), self.Q_star))
        self.C_xpost = numpy.dot((numpy.dot(self.Q, self.C_post)), self.Q.transpose())

        # for i,row in enumerate(self.C_post):
        #    print(list(row))

        # for i, row in enumerate(self.C_xpost):
        #    print(list(row))

        # posterior stds
        self.stds_log_inc  = self.extract_cpost_inc()
        self.stds_log_post = self.extract_cpost()
        
        # posterior mean vector
        if self.pseudo_used:
            self.q_post = numpy.dot(self.C_post, numpy.dot(numpy.dot(self.Q.transpose(), self.C_prior_inv), self.q_prior) + numpy.dot(numpy.dot(Q_star_trans, self.C_x_inv), self.x_star))
        else:
            self.q_post = numpy.dot(self.C_post, (numpy.dot(numpy.dot(Q_star_trans, self.C_x_inv), self.x_star)+numpy.dot(self.C_prior_inv, self.q_prior)))
        self.x_post = numpy.dot(self.Q, self.q_post)
        # print(self.x_post)
        
        # for elem in self.q_post:
        #    print(elem)
        # print('\n')

        # for i,elem in enumerate(self.x_post):
        #    print(elem, ',(', self.quantities[i], ')')

        
    def extract_cpost(self):
        '''
        extract the stds from the posterior diagonal covariance matrix C_post
        '''
        stds = []
        
        for i, row in enumerate(self.C_xpost):
            stds.append(numpy.sqrt(row[i]))

        return stds
    
    def extract_cpost_inc(self):
        '''
        extract the stds from the posterior diagonal covariance matrix C_post
        '''
        stds = []
        for i, row in enumerate(self.C_post):
            stds.append(numpy.sqrt(row[i]))
        return stds

    def build_new_sbtab(self, new_sbtab=False, posterior_sample=False, header=False):
        '''
        generates new SBtab
        '''
        if new_sbtab:
            self.new_rows = []
            for row in new_sbtab: self.new_rows.append(row)
            means  = posterior_sample
            try: self.new_header = header.split('\t')
            except: self.new_header = header
        else:
            means = self.mean_post

        if self.optimized: means = self.mean_post_opt
        else: means = self.mean_post

        finished_rows = [['!!SBtab SBtabVersion="1.0" TableType="Quantity" TableName="Parameter" Document="%s" Date="%s"' % (self.sbtab.filename, datetime.date.today())]]
        finished_rows.append(self.new_header)
        first         = True
        self.hilo     = []
        
        # Add biomass warning if required
        self.check_biomass()
        self.check_enzyme_species()
        
        for i, row in enumerate(self.new_rows):
            if len(row) == len(self.new_header):
                if row[0] == 'QuantityType': continue

                # first: identify the row
                if row[0] in self.species_parameters: row_identifier = (row[0], row[2])
                elif row[0] in self.reaction_parameters: row_identifier = (row[0], row[1])
                else: row_identifier = (row[0], (row[1], row[2]))
                try: row_number = self.id_order[row_identifier]
                except: continue

                # second: fill the row with the balanced values
                if row[0] in self.thermodynamics:
                    row[3] = str(format(float(means[row_number]), '.4f'))
                    row[5] = 'NaN'
                    row[6] = 'NaN'
                else:
                    row[3] = str(format(numpy.exp(float(self.x_post[row_number])), '.4f'))
                    row[5] = str(format(numpy.exp(float(self.x_post[row_number])), '.4f'))
                    row[6] = str(format(numpy.exp(float(self.stds_log_post[row_number])), '.4f'))

                row[7] = str(format(float(means[row_number]), '.4f'))
                row[8] = str(format(float(self.stds_post[row_number] ) , '.4f'))
                if not row in finished_rows:
                    finished_rows.append(row)
                first = self.check_extreme_values(row, first)

        if not first: self.log += '\n'
        if self.hilo != []:
            self.hilo = sorted(self.hilo)
            for entry in self.hilo:
                self.log += entry

                
        return finished_rows

    def check_extreme_values(self, row, first):
        '''
        this function checks whether the given parameter (in form of its corresponding SBtab row)
        has a significantly high or low value.
        '''
        if self.pmin[row[0]] == None and self.pmax[row[0]] == None: return first
        elif float(row[5]) < self.pmin[row[0]]:
            if first:
                self.log += '\n### Warnings about unusually high or low values ### \n'
                first = False
            self.hilo.append('The value for the %s of %s, %s lies under the given lower bound: %s. Please check the accuracy and refer to the FAQ for help.\n'%(row[0], row[1], row[2], row[5]))
        elif float(row[5]) > self.pmax[row[0]]:
            if first:
                self.log += '### Warnings about unusually high or low values ### \n'
                first = False
            self.hilo.append('The value for the %s of %s, %s lies over the given upper bound: %s. Please check the accuracy and refer to the FAQ for help.\n'%(row[0], row[1], row[2], row[5]))

        return first

