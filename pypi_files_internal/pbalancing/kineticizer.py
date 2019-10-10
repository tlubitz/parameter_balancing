#!/usr/bin/env python
import libsbml
import sys
try: from . import misc
except: import misc


class Kineticizer(object):
    '''
    Base kineticizer class. This is a template class. To use it derive it and
    implement the get_denominator function.
    '''
    # mapping of internal names (the ones used in this code) and external names
    # (the ones used in sbtab) the external names have to be changed only at
    # this position and in the sbtab class, the internal names shold *never
    # ever* change
    internal2external = {'km_vec': 'Michaelis constant',
                         'kf': 'substrate catalytic rate constant',
                         'kb': 'product catalytic rate constant',
                         'kv': 'catalytic rate constant geometric mean',
                         'keq': 'equilibrium constant',
                         'met_conc': 'concentration',
                         'enz_conc': 'concentration of enzyme',
                         'mu_vec': 'standard chemical potential',
                         'k_i_vec': 'inhibitory constant',
                         'k_a_vec': 'activation constant',
                         'act_ratio_vec': 'activation baseline ratio',
                         'inh_ratio_vec': 'inhibition baseline ratio',
                         'hill_coeff': 'Hill constant'}

    # mapping of which kinetic parametrisation requires which types of
    # parameters the parameters listed here are automatically retrieved
    # from the sbtab file and written to the sbml (for the parametrisation
    # selected)
    mode2params = {'cat': ('kf', 'kb', 'km_vec', 'k_i_vec', 'k_a_vec',
                           'act_ratio_vec', 'inh_ratio_vec'),
                   'hal': ('kv', 'keq', 'km_vec', 'k_i_vec', 'k_a_vec',
                           'hill_coeff', 'act_ratio_vec', 'inh_ratio_vec'),
                   'weg': ('kv', 'mu_vec', 'km_vec', 'k_i_vec', 'k_a_vec',
                           'hill_coeff', 'act_ratio_vec', 'inh_ratio_vec')}

    # mapping of parameter types and options
    # each parameter option is a tuple of strings
    # 1st: local or global: indicates whether the parameter is associated with
    #      a reaction (local) or not (global)
    # 2nd: '' or 's' or 'm': indicates whether the parameter belongs to a
    #      reaction (local param. and ''), or to a species participating in a
    #      reaction (local param and 's'), to a modifier in a reaction (local
    #      param and 'm'), or just to a species (global param and 's')
    # 3rd: parameter name template ($REAC$ will be replaced by reaction id,
    #      $SPECIES$ will be replaced by species id)
    # 4th: opt or req: indication whether parameter is optional (opt) or
    #      required (req)
    param2options = {'kf': ('local', '', 'kcrf_$REAC$', 'req'),
                     'kb': ('local', '', 'kcrr_$REAC$', 'req'),
                     'km_vec': ('local', 's', 'kmc_$REAC$_$SPECIES$', 'req'),
                     'kv': ('local', '', 'kcrg_$REAC$', 'req'),
                     'keq': ('local', '', 'keq_$REAC$', 'req'),
                     'mu_vec': ('global', 's', 'scp_$SPECIES$', 'req'),
                     'k_i_vec': ('local', 'k', 'kic_$REAC$_$SPECIES$', 'opt'),
                     'k_a_vec': ('local', 'k', 'kac_$REAC$_$SPECIES$', 'opt'),
                     'act_ratio_vec': ('local', 'm', 'rac_$REAC$_$SPECIES$',
                                       'opt'),
                     'inh_ratio_vec': ('local', 'm', 'ric_$REAC$_$SPECIES$',
                                       'opt'),
                     'hill_coeff': ('local', '', 'hco_$REAC$', 'opt'),
                     'met_conc': ('', '', '', 'opt'),
                     'enz_conc': ('', '', '', 'opt')}

    sbo_inh = [20, 206, 207, 536, 537]
    sbo_act = [13, 21, 459, 461, 462]

    type2sbo = {'partial_act': 462,
                'partial_inh': 536,
                'complete_act': 543,
                'complete_inh': 20,
                'specific_act': 533,
                'specific_inh': 206,
                'kf': 320,
                'kb': 321,
                'km_vec': 27,
                'kv': 482,
                'mu_vec': 463,
                'k_i_vec': 261,
                'k_a_vec': 363,
                'act_ratio_vec': 488,
                'inh_ratio_vec': 489,
                'hill_coeff': 382}
    '''
    Class Constructor
    @type  sbtab: libSBtab.Table.KineticDataTable
    @param sbtab: sbtab containg the parameters
    @type  mode: string
    @param mode: type of parametrisation. one of ('cat'|'weg'|'hal')
    @type  enzyme_prefac: boolean
    @param enzyme_prefac: indicate whether enzymes should be included as
                          prefactors to the kinetic laws (true by default)
    @type  default_inh: string
    @param default_inh: default type of inhibition (in case inh. is not
                        specified by SBO). can be ('non-comp', 'comp') for
                        non-competetive / competetive
    @type  default_act: string
    @param default_act: default type of activation. can be ('essential',
                        'non-essntial')
    @type  overwrite_existing: boolean
    @param overwrite_existing: indicate whether existing kinetic laws should
                               be overwritten
    @type  writer: class
    @param writer: writer class providing the write function. default
                   is sys.stderr
    '''
    def __init__(self, model,
                     sbtab=None,
                     mode='cat',
                     enzyme_prefac=True,
                     default_inh='complete_inh',
                     default_act='complete_act',
                     overwrite_existing=True,
                     writer=sys.stderr):

        self._model = model
        self._writer = writer
        if (model.getLevel() < 2) or (model.getLevel() == 2 and
                                      model.getVersion() < 2):
            self._writer.write('SBO Term support requires SBML L2V2 at least. \
SBO Terms not supported for this document.\n')
        self._sbtab = sbtab
        self._enzyme_prefac = enzyme_prefac

        if default_inh not in ['partial_inh', 'complete_inh', 'specific_inh']:
            raise Exception('Unknown inhibition')
        self._default_inh = default_inh

        if default_act not in ['partial_act', 'complete_act', 'specific_act']:
            raise Exception('Unknown activation')
        self._default_act = default_act

        self.sbo2type = dict(zip(self.type2sbo.values(), self.type2sbo.keys()))

        if sbtab is not None:
            [global_params, local_params] = self._pack_parameters(mode)
        else:
            [global_params, local_params] = self._default_parameters(mode)

        # write global parameters to SBML
        self._create_global_params(global_params, mode)

        # write the kinetic laws
        for i, r in enumerate(model.getListOfReactions()):
            if r.getKineticLaw() is not None and not overwrite_existing:
                continue
            self._assign_kinetic(r, local_params[i], mode)
        self._set_metabolite_concentrations()
        try:
            self._set_enzyme_concentrations()
        except:
            pass

    def _print_warning(self, text):
        '''
        internal function to print a warning to the standard output
        (writer class)
        '''
        self._writer.write(text + '\n')

    def _assign_kinetic(self, reaction, params, mode):
        '''
        create a kinetic law for the given reaction
        @type  reaction:  libsbml.reaction
        @param reaction:  the reaction
        @type  params:    list
        @param params:    list of local (assigned to a reaction) parameters,
                          as generated by pack_parameters function
        @type  mode:      string
        @param mode:      parametrisation type ('cat' | 'weg' | 'hal')
        '''
        kl = reaction.createKineticLaw()
        kl.setSBOTerm(self._kinetic_law_sbo)
        # create the local parameters
        self._create_local_params(reaction, params, mode)

        # generate numerator, denominator, regulational prefactors
        numer = self._get_numerator(reaction, mode)
        denom = self._get_denominator(reaction, mode)
        prefac, denom_term = self._get_regulation(reaction)

        # add regulation to denominator, if there is
        if denom_term != '':
            denom = '(%s + %s)' % (denom_term, denom)

        # generate the kinetic formula
        formula = '%s / %s' % (numer, denom)
        if prefac != '':
            formula = '(%s) * %s' % (prefac, formula)

        if self._enzyme_prefac:
            # include enzyme as prefactor to kinetic law
            concentration = self._get_sbtab_entry('enz_conc', reaction)
            compartment = 1.
            enz_fac = '( %s / %s)' % (concentration, compartment)
            formula = enz_fac + ' * ' + formula

        kl.setFormula(formula)

    def _add_param(self, p_type, value, reaction=None, species=None):
        '''
        generic function to add a parameter to a reaction (local param.) of to
        the sbml model (global param)
        @type  p_type:     string
        @param p_type:     type of the parameter (key of self.param2options)
        @type  value:      float
        @param value:      value
        @type  reaction:   libsbml.reaction
        @param reaction:   the reaction, if paramter is associated with
                           a reaction
        @type  species:    libsbml.species
        @param species:    the spiecies if parameter is associated with
                           a species
        '''
        # get the options for this parameter
        scope, deps, id_temp, required = self.param2options[p_type]
        # add reaction_id and species_id to parameter id, if possible
        if reaction:
            id_temp = id_temp.replace('$REAC$', reaction.getId())
        if species:
            try: s_id = species.getSpecies()
            except: s_id = species.getId()
            id_temp = id_temp.replace('$SPECIES$', s_id)

        # add parameter to reaction (local scope) of to mode (global scope)
        if scope == 'local':
            kl = reaction.getKineticLaw()
            # before creating a parameter,
            # check if such a parameter already exists
            if id_temp in [p.getId() for p in kl.getListOfParameters()]:
                return
            p = kl.createParameter()
        else:
            # before creating a parameter,
            # check if such a parameter already exists
            if id_temp in [p.getId() for p in self._model.getListOfParameters()]:
                return
            p = self._model.createParameter()

        p.setId(id_temp)
        p.setName(id_temp)

        try:
            p.setSBOTerm(self.type2sbo[p_type])
        except: pass

        p.setValue(float(value))

    def _create_local_params(self, reaction, params, mode):
        '''
        create the local parameters for a reaction
        @type  reaction:   libsbml.reaction
        @param reaction:   the reaction, if paramter is associated with
                           a reaction
        @type  params:     list
        @param params:     list of local (assigned to a reaction) parameters,
                           as generated by pack_parameters function
        @type  mode:       string
        @param mode:       parametrisation type ('cat' | 'weg' | 'hal')
        '''
        local_param_count = 0

        # iterate over all parameter types for this mode
        for p_type in self.mode2params[mode]:
            # get options
            scope, deps, id_temp, required = self.param2options[p_type]
            if scope == 'global':  # this function is for global parameters
                continue
            # get the value
            p_value = params[local_param_count]
            local_param_count += 1
            # based on dependencies,
            # add one parameter for reaction / reactant / modifier
            if deps == '':
                self._add_param(p_type, p_value, reaction)
            elif deps == 's':
                for i, s in enumerate(misc.get_participants(reaction)):
                    self._add_param(p_type, p_value[i], reaction, s)
            elif deps == 'k':      # here go the k_i_vec and k_a_vec
                for i, s in enumerate(misc.get_modifiers(reaction)):
                    sboterm = s.getSBOTerm()
                    if p_value[i] is None:
                        continue
                    else:
                        if sboterm in self.sbo_act and p_type == 'k_a_vec':
                            self._add_param(p_type, p_value[i], reaction, s)
                        elif sboterm in self.sbo_inh and p_type == 'k_i_vec':
                            self._add_param(p_type, p_value[i], reaction, s)
            elif deps == 'm':     # here go the inh_ratio_vec and act_ratio_vec
                for i, s in enumerate(reaction.getListOfModifiers()):
                    if p_value[i] is None:
                        continue
                    self._add_param(p_type, p_value[i], reaction=reaction,
                                    species=s)

    def _create_global_params(self, g_params, mode):
        '''
        create the global parameters
        @type  g_params:  list
        @param g_params:  list of global parameters,  as generated by
                          pack_parameters function
        @type  mode:      string
        @param mode:      parametrisation type ('cat' | 'weg' | 'hal')
        '''
        # Hack, add temperature and gas constant if parametrization is 'weg'
        if mode == 'weg':
            p = self._model.createParameter()
            p.setId('temp')
            p.setValue(300)
            p = self._model.createParameter()
            p.setId('R')
            p.setValue(8.31)

        # iterate over all global parameter types
        global_param_count = 0
        for p_type in self.mode2params[mode]:
            scope, deps, id_temp, required = self.param2options[p_type]
            if scope == 'local':
                continue
            p_value = g_params[global_param_count]
            global_param_count += 1
            # decide whether parameter depends on species / reaction and add it
            if deps == 's' or deps == 'k':
                i = 0
                for s in self._model.getListOfSpecies():
                    if misc.is_enzyme(s):
                        continue
                    self._add_param(p_type, p_value[i], reaction=None,
                                    species=s)
                    i += 1
            elif deps == 'r':
                for i, r in enumerate(self._model.getListOfReactions):
                    self._add_param(p_type, p_value[i], reaction=r,
                                    species=None)

    def _get_sbtab_entry(self, param_type, reaction=None, species=None):
        '''
        get sbtab entry for libsbml reaction and species
        @type  param_type:  string
        @param param_type:  type of parameters (the left side
                            of self.internal2external)
        @type  reaction:    libsbml.reaction
        @param reaction:    reaction
        @type  species:     libsbml.species
        @param species:     species
        @rtype:             float or None
        @return:            return value of None if not found
        '''
        def get_sbtab_entry(sbtab, **arg_dict):
            # generic function to get any sbtab entry
            # iterates over all rows in the sbtab and checks wehter all given
            # keys are right            
            for row in sbtab.value_rows:
                if len(row) == len(sbtab.columns):
                    r_column = ''
                    c_column = ''

                    for key in arg_dict:
                        if key == 'Reaction':
                            r_column = sbtab.columns_dict['!Reaction:SBML:reaction:id']
                        if key == 'QuantityType':
                            qt_column = sbtab.columns_dict['!QuantityType']
                        if key == 'Compound':
                            c_column = sbtab.columns_dict['!Compound:SBML:species:id']
                        m_column = sbtab.columns_dict['!Mode']

                        
                    if len(arg_dict) > 2 and row[qt_column] == arg_dict['QuantityType']:
                        if row[r_column] == arg_dict['Reaction']:
                            if row[c_column] == arg_dict['Compound']:
                                return row[m_column]
                            '''
                            # this here turned out to mess up reactions with multiple modifiers
                            elif row[qt_column] == 'inhibitory constant' or row[qt_column] == 'activation constant':
                                return row[m_column]
                            '''
                        elif arg_dict['Reaction'] is None:
                            if row[c_column] == arg_dict['Compound']:
                                return row[m_column]
                    else:
                        if row[qt_column] == arg_dict['QuantityType'] and row[r_column] == arg_dict['Reaction']:
                            return row[m_column]

        # get the reaction and species IDs
        r_id = s_id = None
        try: r_id = reaction.getId()
        except: pass
        try: s_id = species.getId()
        except: pass

        # check whether this parameter is required
        required = self.param2options[param_type][-1] == 'req'
        
        # get the parameter value
        if param_type == 'hill_coeff': return 1
        elif param_type == 'act_ratio_vec' or param_type == 'inh_ratio_vec':
            return None
        try:
            if s_id is not None:
                try:
                    value = float(get_sbtab_entry(self._sbtab,
                                                  QuantityType=self.internal2external[param_type],
                                                  Reaction=r_id,
                                                  Compound=s_id))
                except:
                    value = None
            else:
                value = float(get_sbtab_entry(self._sbtab,
                                              QuantityType=self.internal2external[param_type],
                                              Reaction=r_id))
        except AttributeError:
            # if no value is given, raise error or continue
            if required:
                raise Exception('''No value found for type %s: (species %s and
                reaction %s)''' % (self.internal2external[param_type],
                                   s_id, r_id))
            value = None

        return value

    def _pack_parameters(self, mode):
        '''
        get the parameters from the sbtab and pack them in a way that they can
        be added easily. Parameters are separated in local ones (associated
        with a reaction) and global ones.
        @type  mode:  string
        @param mode:  parametrisation type ('cat' | 'weg' | 'hal')
        @rtype:       tuple
        @return:      tuple of (global_parameters, local_parameters) where
                      global parameters is the list of global params (1 entry
                      for each param type) and local_parameter is the list of
                      local params (1 entry for each reaction)
        '''
        # get the global params
        global_params = []
        # get global params for all parameters specified for this mode
        for p_type in self.mode2params[mode]:
            p_vec = []
            scope, deps, dummy, required = self.param2options[p_type]
            if scope == 'local':
                continue
            # if param. depends on species, get one param for each species, if
            # it depends on reaction get one for each reaction
            if deps == 's' or deps == 'k':
                for s in self._model.getListOfSpecies():
                    if misc.is_enzyme(s):
                        continue
                    p_value = self._get_sbtab_entry(p_type, species=s)
                    p_vec.append(p_value)
            elif deps == 'r':
                for r in self._model.getListOfReactions():
                    p_value = self._get_sbtab_entry(p_type, reaction=r)
                    p_vec.append(p_value)
            global_params.append(p_vec)

        # get local parameters (one set for each reaction)
        local_params = []
        for r in self._model.getListOfReactions():
            reaction_params = []
            for p_type in self.mode2params[mode]:
                scope, deps, dummy, required = self.param2options[p_type]
                if scope == 'global':
                    continue
                if deps == '':
                    p_vec = self._get_sbtab_entry(p_type, reaction=r)
                elif deps == 's':
                    p_vec = []
                    for s in misc.get_participants(r):
                        p_value = self._get_sbtab_entry(p_type, reaction=r,
                                                        species=self._model.getSpecies(s.getSpecies()))
                        p_vec.append(p_value)
                elif deps == 'k':
                    p_vec = []
                    for s in misc.get_modifiers(r):

                                            
                        p_value = self._get_sbtab_entry(p_type, reaction=r,
                                                        species=self._model.getSpecies(s.getSpecies()))
                        p_vec.append(p_value)
                elif deps == 'm':
                    p_vec = []
                    for m in r.getListOfModifiers():
                        p_value = self._get_sbtab_entry(p_type, reaction=r,
                                                        species=self._model.getSpecies(m.getSpecies()))
                        p_vec.append(p_value)
                reaction_params.append(p_vec)
            local_params.append(reaction_params)
            
        return (global_params, local_params)

    def _default_parameters(self, mode):
        '''
        if not sbtab is given, this function provides default parameters
        @type  mode:  string
        @param mode:  parameparametrisation type ('cat' | 'weg' | 'hal')
        @rtype:       tuple
        @return:      tuple of (global_parameters,local_parameters) where
                      global parameters is the list of global params (1
                      entry for each param type) and local_parameter is
                      the list of local params (1 entry for each reaction)
        '''
        self._writer.write('No Parameters given. Setting all to 1.')
        global_params = []
        local_params = []
        # generate all required local parameters
        for p_type in self.mode2params[mode]:
            scope, deps, dummy, required = self.param2options[p_type]
            if scope != 'global':
                continue
            if deps == 's' or deps == 'k':
                global_params.append([1.] * self._model.getNumSpecies())
            elif deps == 'r':
                global_params.append([1.] * self._model.getNumReactions())

        # generate all required global parameters
        for r in self._model.getListOfReactions():
            if scope != 'local':
                continue
            reaction_params = []
            for p_type in self.mode2params[mode]:
                scope, deps, dummy, required = self.param2options[p_type]
                if deps == '':
                    reaction_params.append(1.)
                elif deps == 's':
                    reaction_params.append([1.] * (r.getNumReactants() +
                                                   r.getNumProducts()))
                elif deps == 'k':
                    reaction_params.append([1] * r.getNumModifiers())
                elif deps == 'm':
                    reaction_params.append([1] * r.getNumModifiers())
            local_params.append(reaction_params)
        return (global_params, local_params)

    def _set_metabolite_concentrations(self):
        '''
        set the metabolite concentrations to values given in the sbtab
        (if specified)
        '''
        for s in self._model.getListOfSpecies():
            if misc.is_enzyme(s):  # enzymes are treated separately
                continue
            conc = float(self._get_sbtab_entry('met_conc', species=s))
            try:
                conc = float(self._get_sbtab_entry('met_conc', species=s))
                s.setInitialConcentration(conc)
            except:
                self._print_warning('''No initial concentration found for
                species %s. Setting value to 1.''' % s.getId())
                s.setInitialConcentration(1.)

    def _set_enzyme_concentrations(self):
        '''
        set the enzyme concentrations to values given in the sbtab
        '''
        for r in self._model.getListOfReactions():
            enz = misc.get_enzyme_for_reaction(r)
            try:
                conc = float(self._get_sbtab_entry('enz_conc', reaction=r))
                enz.setInitialConcentration(conc)
            except:
                self._print_warning('''No concentration found for enzyme %s.
                Setting value to 1.''' % enz.getId())
                enz.setInitialConcentration(1.)

    def _get_substrate_term(self, species_ref):
        '''
        get term for species as it appears in the kinetic law
        (concentration/KM)^stoich
        @type  species_ref:  libsbml.speciesReference
        @param species_ref:  species reference from reaction
        @rtype:              string
        @return:             (conc./KM)^stoich.
        '''
        reaction = species_ref.getParentSBMLObject().getParentSBMLObject()
        km = self.param2options['km_vec'][2].replace('$REAC$',
                                                     reaction.getId()).replace('$SPECIES$',
                                                                               species_ref.getSpecies())
        try: stoich = species_ref.getStoichiometry()
        except: stoich = '1'
        text = '(%s/%s)' % (species_ref.getSpecies(), km)
        if stoich != 1:
            text += '^' + str(stoich)
        return text

    def _get_numerator(self, reaction, mode):
        '''
        get the numerator for the kinetic law
        @type  reaction:  libsbml.reaction
        @param reaction:  the reaction
        @type  mode:      string
        @param mode:      parametrisation type ('cat' | 'weg' | 'hal')
        '''
        # define little functionto get s^stoich
        def get_substrate_term(s):
            s_id = s.getSpecies()
            if s.getStoichiometry() != 1:
                s_id += '^' + str(s.getStoichiometry())
            return s_id

        # define little function to get KM^stoich
        def get_km_term(s):
            km_id = self.param2options['km_vec'][2].replace('$REAC$',
                                                            reaction.getId()).replace('$SPECIES$',
                                                                                      s.getSpecies())
            if s.getStoichiometry() != 1:
                km_id += '^' + str(s.getStoichiometry())
            return km_id

        # get the id templates for kf, kb and hill coeff.
        kf_id_temp = self.param2options['kf'][2]
        kb_id_temp = self.param2options['kb'][2]
        hill_id_temp = self.param2options['hill_coeff'][2]

        # generate the numerator term dpending on parameterization type
        if mode == 'cat':
            forw_term = kf_id_temp.replace('$REAC$', reaction.getId()) + ' * '
            back_term = kb_id_temp.replace('$REAC$', reaction.getId()) + ' * '
            forw_term += ' * '.join([self._get_substrate_term(s)
                                     for s in reaction.getListOfReactants()])
            back_term += ' * '.join([self._get_substrate_term(s)
                                     for s in reaction.getListOfProducts()])
            if reaction.getNumReactants() == 0:
                forw_term = '0'
            elif reaction.getNumProducts() == 0:
                back_term = '0'
            term = '( ( %s ) - ( %s ) )' % (forw_term, back_term)
        elif mode == 'hal':
            react_term = ' * '.join([get_substrate_term(s)
                                      for s in reaction.getListOfReactants()]) \
                                          or '0'
            prod_term = ' * '.join([get_substrate_term(s)
                                     for s in reaction.getListOfProducts()]) \
                                         or '0'
            km_term = ' * '.join([get_km_term(s)
                                    for s in misc.get_participants(reaction)])
            kv = self.param2options['kv'][2].replace('$REAC$',
                                                     reaction.getId())
            keq = self.param2options['keq'][2].replace('$REAC$',
                                                       reaction.getId())
            hill_coeff = hill_id_temp.replace('$REAC$', reaction.getId())
            term = '''( %s * ( ( ((%s)^(%s/2.)) * %s) - (((%s)^(-%s/2.)) * %s)
            ) )''' % (kv, keq, hill_coeff, react_term, keq, hill_coeff,
                      prod_term)
            term = '(%s) / sqrt(%s)' % (term, km_term)
        elif mode == 'weg':
            mu_id_temp = self.param2options['mu_vec'][2]

            def get_mu_term(s):
                mu_id = mu_id_temp.replace('$SPECIES$', s.getSpecies())
                if s.getStoichiometry() != 1:
                    mu_id = str(s.getStoichiometry()) + '*' + mu_id
                return mu_id

            mu_f_term = ' + '.join([get_mu_term(s)
                                    for s in reaction.getListOfReactants()]) \
                                        or '0'
            mu_b_term = ' + '.join([get_mu_term(s)
                                    for s in reaction.getListOfProducts()]) \
                                        or '0'
            mu_term = '(%s) - (%s)' % (mu_b_term, mu_f_term)
            hill_coeff = hill_id_temp.replace('$REAC$', reaction.getId())
            react_term = ' * '.join([get_substrate_term(s)
                                      for s in reaction.getListOfReactants()])\
                                          or '0'
            prod_term = ' * '.join([get_substrate_term(s)
                                     for s in reaction.getListOfProducts()]) \
                                         or '0'
            km_term = ' * '.join([get_km_term(s)
                                   for s in misc.get_participants(reaction)])
            kv = self.param2options['kv'][2].replace('$REAC$',
                                                     reaction.getId())
            term = '''( %s * (( exp(-%s * (%s)/(2*R*temp)) * (%s) ) -
            (exp(%s * (%s)/(2*R*temp)) * (%s)) ) / sqrt(%s))''' \
                % (kv, hill_coeff, mu_term, react_term, hill_coeff, mu_term,
                   prod_term, km_term)
        else:
            # in case somebody chose a mode we do not know
            raise Exception('Implement')

        return term

    def _get_denominator(self, reaction, mode):
        '''
        this class is just a template. this function needs to be implemented in
        the derived class
        @type  reaction: libsbml.reacion
        @param reaction: the reaction
        @type  mode:     string
        @param mode:     parameparametrisation type ('cat' | 'weg' | 'hal')
        @rtype:          string
        @return:         formula for the denominator
        '''
        raise Exception('Implement me')

    def _get_regulation(self, reaction):
        '''
        get regulational term for a reaction
        @type  reaction:  libsbml.reaction
        @param reaction:  reaction
        @rtype:           tuple
        @return:          (prefactor, denominator_term), where prefactor is a
                          list of prefactors to the kinetic law and denominator
                          term is a list of terms added to the denominator (one
                          list entry per modifier)
        '''
        prefacs = []
        denom_terms = []
        # generate one term for each modifier
        for m in reaction.getListOfModifiers():
            if misc.is_enzyme(self._model.getSpecies(m.getSpecies())):
                continue
            # get the ids for the possiby involved parameters
            ka = self.param2options['k_a_vec'][2].replace('$REAC$',
                                                          reaction.getId()).replace('$SPECIES$',
                                                                                    m.getSpecies())
            ki = self.param2options['k_i_vec'][2].replace('$REAC$',
                                                          reaction.getId()).replace('$SPECIES$',
                                                                                    m.getSpecies())
            act_ratio = self.param2options['act_ratio_vec'][2].replace('$REAC$',
                                                                       reaction.getId()).replace('$SPECIES$',
                                                                                                 m.getSpecies())
            inh_ratio = self.param2options['inh_ratio_vec'][2].replace('$REAC$',
                                                                       reaction.getId()).replace('$SPECIES$',
                                                                                                 m.getSpecies())
            act_term = '(%s/%s)' % (m.getSpecies(), ka)
            inh_term = '(%s/%s)' % (m.getSpecies(), ki)

            # decide for regulation type
            sbo = m.getSBOTerm()
            try:
                reg_type = self.sbo2type[sbo]
            except:
                if reaction.getKineticLaw().getParameter(ka):
                    reg_type = self._default_act
                elif reaction.getKineticLaw().getParameter(ki):
                    reg_type = self._default_inh
                else:
                    continue

            if 'act' in reg_type:
                if not reaction.getKineticLaw().getParameter(ka):
                    raise Exception('Parameter %s not given' % ka)
            if 'inh' in reg_type:
                if not reaction.getKineticLaw().getParameter(ki):
                    raise Exception('Parameter %s not given' % ki)
            if reg_type == 'partial_act' \
               and not reaction.getKineticLaw().getParameter(act_ratio):
                    raise Exception('Parameter %s not given' % act_ratio)
            if reg_type == 'partial_inh' \
               and not reaction.getKineticLaw().getParameter(inh_ratio):
                    raise Exception('Parameter %s not given' % inh_ratio)
            reaction.getKineticLaw().setSBOTerm(self.type2sbo[reg_type])

            # generate the actual terms
            if reg_type == 'partial_act':
                term = '(%s + (1 - %s) * ( (%s) / (1 + %s)))' % (act_ratio,
                                                                 act_ratio,
                                                                 act_term,
                                                                 act_term)
            elif reg_type == 'partial_inh':
                term = '(%s + (1 - %s) * ( 1 / (1 + %s)))' % (inh_ratio,
                                                              inh_ratio,
                                                              inh_term)
            elif reg_type == 'complete_act':
                term = '((%s)/(1 + %s))' % (act_term, act_term)
            elif reg_type == 'complete_inh':
                term = '(1/(1 + %s))' % (inh_term)
            elif reg_type == 'specific_act':
                t = '(%s/%s)' % (ka, m.getSpecies())
                denom_terms.append(t)
                continue
            elif reg_type == 'specific_inh':
                denom_terms.append(inh_term)
                continue
            else:
                raise Exception("Unknown regulation type: typo in code")
            prefacs.append(term)

        prefac = ' * '.join(prefacs)
        denom_term = ' + '.join(denom_terms)

        return (prefac, denom_term)


class KineticizerCS(Kineticizer):
    ''' common saturable implementation of kineticizer '''
    _kinetic_law_sbo = 528

    def _get_denominator(self, reaction, mode):
        get_substr_term = lambda s: '(' + self._get_substrate_term(s).replace('(', '1 + (') + ')'
        forw = ' * '.join([get_substr_term(s)
                           for s in reaction.getListOfReactants()])
        backw = ' * '.join([get_substr_term(s)
                            for s in reaction.getListOfProducts()])
        if not backw:
            backw = '0'
        if not forw:
            forw = '0'
        text = '(( %s ) + ( %s ) - 1)' % (forw, backw)

        return text


class KineticizerMS(Kineticizer):
    ''' multiplicative saturable implementation of kineticizer'''
    _kinetic_law_sbo = 530

    def _get_denominator(self, reaction, mode):
        get_substr_term = lambda s: self._get_substrate_term(s).replace('(', '( 1 + ')
        text = ' * '.join([get_substr_term(s)
                           for s in misc.get_participants(reaction)])
        return '(%s)' % text


class KineticizerDS(Kineticizer):
    ''' direct saturable implementation '''
    _kinetic_law_sbo = 529

    def _get_denominator(self, reaction, mode):
        forw = ' * '.join([self._get_substrate_term(s)
                           for s in reaction.getListOfReactants()])
        backw = ' * '.join([self._get_substrate_term(s)
                            for s in reaction.getListOfProducts()])
        if not forw:
            forw = '0'
        if not backw:
            backw = '0'
        text = '( 1 + %s + %s )' % (forw, backw)
        return text


class KineticizerFD(Kineticizer):
    ''' fd implementation '''
    _kinetic_law_sbo = 532

    def _get_denominator(self, reaction, mode):

        def get_substr_term(s):
            term = self._get_substrate_term(s)
            term = term[:term.rfind(')') + 1]
            term = 'sqrt(%s)' % term
            return term

        text = ' * '.join([get_substr_term(s)
                           for s in misc.get_participants(reaction)])
        return '(%s)' % text


class KineticizerRP(Kineticizer):
    ''' mass action implementation '''
    _kinetic_law_sbo = 531


    def _get_denominator(self, reaction, mode):
        return '(1)'


if __name__ == '__main__':

    d = libsbml.readSBML(sys.argv[1])
    m = d.getModel()

    sbtab = libSBtab.Table.KineticDataTable()
    sbtab.fromTSVFile(sys.argv[2])

    KineticizerRP(m, sbtab, mode='cat')
    print('<?xml version="1.0" encoding="UTF-8"?>\n' + d.toSBML())
