# -*- coding: utf-8 -*-
# this file is released under public domain and you can use without limitations
#########################################################################
## - download is for downloading files uploaded in the db (does streaming)
## - call exposes all registered services (none by default)
#########################################################################

# some packages that need to be imported
# libsbml needs to be installed, the rest are part of the repository
import libsbml
import misc
import copy
import balancer
import kineticizer
import SBtab
import tablibIO
import validatorSBtab
import sbtab2sbml
import sbml2sbtab


def index():
    '''
    index function that is called upon first call of website
    '''
    redirect(URL('../../static/css/css_template/main.html'))


def clearsession():
    '''
    function to clear all session variables in case that anything goes awry in
    the session. It basically just resets all session variables to empty,
    False, or None, depending on the type.
    '''
    session.warnings_sbml = []
    session.warnings_sbtab = []
    session.warnings_prior = []
    session.warnings_fl = []
    session.sbml = None
    session.sbml_name = None
    session.result_sbml = None
    session.sbmls = []
    session.sbml_names = []
    session.result_sbml_name = None
    session.sbtab = None
    session.sbtab_name = None
    session.sbtab_fls = []
    session.sbtab_fl_names = []
    session.sbtab_fl = None
    session.sbtab_fl_name = None
    session.sbtabs = []
    session.result_sbtab = []
    session.result_sbtab_name = []
    session.sbtab_names = []
    session.prior = None
    session.prior_name = None
    session.priors = []
    session.prior_names = []
    session.emptysbtab = False
    session.log = None
    session.config_file = None
    session.config_filename = None
    session.parameter_dict = {}
    redirect('http://www.parameterbalancing.net')


def balancing():
    '''
    first page for the online balancing (see def balancing_old for
    intermediate interface): upload of files
    '''
    response.title = T('Parameter Balancing for Kinetic Models of Cell '
                       'Metabolism')
    response.subtitle = T('Online Balancing')

    ##########################################################################
    ### LEFT SIDE: SBML (mandatory), SBtab (opt), Prior (opt), Config (opt)###
    ##########################################################################
    ###                          1: SBML                                  ####
    ##########################################################################

    sbmlform = SQLFORM.factory(Field('File', 'upload', uploadfolder="/tmp",
                                     label='Upload SBML file (.xml)',
                                     requires=IS_LENGTH(10242880, 1,
                                                        error_message='Max upload size: 5MB')))

    # file is uploaded, checked, and added to the list; variables are declared
    if sbmlform.process(formname='form_one').accepted:
        response.flash = 'form accepted'
        session.warnings_sbml = []
        if 'sbmls' not in session:
            session.sbmls = []
            session.sbml_names = []
        filename = request.vars.File.filename
        if filename in session.sbml_names:
            session.warnings_sbml.append('Error: Duplicate file name. Please '
                                         'remove the file with the same name '
                                         'before uploading.')

        try:
            valid_extension = misc.validate_file_extension(filename, 'sbml')
            if not valid_extension:
                session.warnings_sbml.append('Error: The file extension of '
                                             'file %s is not compliant with '
                                             'SBML. Please upload an SBML '
                                             'model with the .xml extension'
                                             '.' % (filename))
            else:
                session.sbml = request.vars.File.value.decode('utf-8')
                session.sbml_name = filename
                session.sbmls.append(session.sbml)
                session.sbml_names.append(filename)
                size_warning = misc.size_warning(session.sbml)
                if size_warning:
                    session.warnings_sbml.append(size_warning)
        except:
            session.warnings_sbml.append('Error: This is no valid SBML file')
        try: redirect(URL('../../default/balancing'))
        except: redirect(URL('../balancing'))

    elif sbmlform.errors: response.flash = 'Form has errors'

    ###########################################################################
    ###                          2: SBtab                                 #####
    ###########################################################################

    sbtabform = SQLFORM.factory(Field('File', 'upload',
                                      uploadfolder="/tmp",
                                      label='Upload data file (.tsv)',
                                      requires=IS_LENGTH(10485760, 1,
                                                         error_message='Max upload size: 10MB')))
    # upload form
    if sbtabform.process(formname='form_two').accepted:
        response.flash = 'form accepted'
        session.warnings_sbtab = []
        if 'sbtabs' not in session:
            session.sbtabs = []
            session.sbtab_names = []
        filename = request.vars.File.filename
        if filename in session.sbtab_names:
            session.warnings_sbtab.append('Error: Duplicate file name. Please '
                                          'remove the file with the same name '
                                          'before uploading.')
        sbtab_file = request.vars.File.value
        try:
            sbtab_data = SBtab.SBtabTable(sbtab_file.decode('utf-8'), filename)
            if sbtab_data.table_type != 'Quantity':
                session.warnings_sbtab.append('Error: The data file has an inc'
                                              'orrect table type: "'"%s"'". Th'
                                              'e correct table type would '
                                              'be "'"Quantity"'".\n' % (sbtab_data.table_type))
            valid_extension = misc.validate_file_extension(request.vars.File.filename,
                                                           'sbtab')
            if not valid_extension:
                session.warnings_sbtab.append('Error: The file extension is no'
                                              't compliant with SBtab. Please '
                                              'upload an SBtab parameter file '
                                              'with the .tsv extension.')
            else:
                try:
                    session.sbtab = sbtab_data
                    session.sbtab_name = sbtab_data.filename
                    session.sbtabs.append(sbtab_data)
                    session.sbtab_names.append(sbtab_data.filename)
                    # upon upload we check whether the SBtab file holds any IDs
                    # which cannot be found in the SBML file
                    if session.sbml:
                        reader = libsbml.SBMLReader()
                        sbml = reader.readSBMLFromString(session.sbml)
                        sbml_model = sbml.getModel()
                        sbtabid2sbmlid = misc.id_checker(sbtab_data,
                                                         sbml_model)
                        if sbtabid2sbmlid != []:
                            for element in sbtabid2sbmlid:
                                session.warnings_sbtab.append(element)
                except:
                    session.warnings_sbtab.append('Error: The SBtab file does '
                                                  'not fit to the SBML file.')
        except:
            session.warnings_sbtab.append('This is not a valid SBtab file. Ple'
                                          'ase ensure the validity of the data'
                                          ' file.')

        try:
            def_file_open = open('./applications/pb/static/files/default_files'
                                 '/definitions.tsv')
            definition_string = def_file_open.read()
            definition_name = 'definitions.tsv'
            sbtab_def = SBtab.SBtabTable(definition_string, definition_name)
            TableValidClass = validatorSBtab.ValidateTable(sbtab_data,
                                                           sbtab_def)
            warnings = TableValidClass.return_output()
            if warnings != []:
                for w in warnings:
                    session.warnings_sbtab.append(w)
        except:
            session.warnings_sbtab.append('Error: The data file could not be v'
                                          'alidated. We suggest to validate th'
                                          'e file manually on www.sbtab.net.')
    elif sbtabform.errors: response.flash = 'Form has errors'

    ###########################################################################
    ###                          3: PRIOR                                 #####
    ###########################################################################

    priorform = SQLFORM.factory(Field('File', 'upload', uploadfolder="/tmp",
                                      label='Upload prior file (.tsv)',
                                      requires=IS_LENGTH(10485760, 1,
                                                         error_message='Max upload size: 10MB')))

    # upload form
    if priorform.process(formname='form_three').accepted:
        response.flash = 'form accepted'
        session.warnings_prior = []
        if 'prior_names' not in session:
            session.priors = []
            session.prior_names = []
        filename = request.vars.File.filename
        if filename in session.prior_names:
            session.warnings_prior.append('Error: Duplicate file name. Please '
                                          'remove the file with the same name '
                                          'before uploading.')

        try:
            prior_file = request.vars.File.value
            sbtab_prior = SBtab.SBtabTable(prior_file.decode('utf-8'),
                                           filename)
            if sbtab_prior.table_type != 'Quantity':
                session.warnings_prior.append('Error: The data file has an inc'
                                              'orrect table type: "'"%s"'". Th'
                                              'e correct table type would be'
                                              '"'"Quantity"'".\n' % (sbtab_prior.table_type))
        except:
            session.warnings_prior.append('Error: This is not a valid SBtab pr'
                                          'ior file.')

        try:
            def_file_open = open('./applications/pb/static/files/default_files'
                                 '/definitions.tsv')
            definition_file = def_file_open.read()
            definition_name = 'definitions.tsv'
            sbtab_def = SBtab.SBtabTable(definition_file, definition_name)
            TableValidClass = validatorSBtab.ValidateTable(sbtab_prior,
                                                           sbtab_def)
            warnings = TableValidClass.return_output()
            if warnings != []:
                session.warnings_prior.append(w)
        except:
            session.warnings_prior.append('Error: The prior data file could no'
                                          't be validated. We suggest to valid'
                                          'ate the file manually on www.sbtab.'
                                          'net.')

        # check file extension
        try:
            valid_extension = misc.validate_file_extension(sbtab_prior.filename,
                                                           'sbtab')
            if not valid_extension:
                session.warnings_prior.append('Error: The file extension is no'
                                              't compliant with the prior tabl'
                                              'e. Please upload a file with th'
                                              'e .tsv extension.')
            else:
                try:
                    prior = sbtab_prior
                    delimiter = misc.check_delimiter(sbtab_prior.return_table_string())
                    validity = misc.valid_prior(sbtab_prior)
                    for warning in validity:
                        session.warnings_prior.append(warning)
                    session.priors.append(sbtab_prior)
                    session.prior_names.append(sbtab_prior.filename)
                except:
                    session.warnings_prior.append('Error: The prior file is no'
                                                  't valid')
        except:
            session.warnings_prior.append('Error: This is not a valid SBtab pr'
                                          'ior file.')
    elif priorform.errors: response.flash = 'Form has errors'

    ###########################################################################
    ###                          4: CONFIG                                #####
    ###########################################################################

    configform = SQLFORM.factory(Field('File', 'upload', uploadfolder="/tmp",
                                       label='Upload options file (.tsv)',
                                       requires=IS_LENGTH(10485760, 1,
                                                          error_message='Max upload size: 10MB')))

    # upload form
    if configform.process(formname='form_four').accepted:
        response.flash = 'form accepted'
        session.warnings_config = []
        filename = request.vars.File.filename
        config_file = request.vars.File.value
        try:
            sbtab_config = SBtab.SBtabTable(config_file.decode('utf-8'),
                                            filename)
            valid_extension = misc.validate_file_extension(sbtab_config.filename,
                                                           'sbtab')
            if sbtab_config.table_type != 'PbConfig':
                session.warnings_config.append('Error: The SBtab file has an i'
                                               'ncorrect table type: "'"%s"'".'
                                               ' The correct table type would '
                                               'be "'"PbConfig"'".\n' % (TableValidClass.sbtab.table_type))
            if not valid_extension:
                session.warnings_config.append('Error: The file extension of f'
                                               'ile %s is not compliant with S'
                                               'Btab. Please upload a options '
                                               'SBtab file with the .tsv exten'
                                               'sion.' % (filename))
            else:
                session.config_file = sbtab_config
                session.config_filename = sbtab_config.filename
                try:
                    def_file_open = open('./applications/pb/static/files/defau'
                                         'lt_files/definitions.tsv')
                    definition_file = def_file_open.read()
                    definition_name = 'definitions.tsv'
                    sbtab_def = SBtab.SBtabTable(definition_file,
                                                 definition_name)
                    TableValidClass = validatorSBtab.ValidateTable(sbtab_config,
                                                                   sbtab_def)
                    warnings = TableValidClass.return_output()
                    if warnings != []:
                        for w in warnings:
                            session.warnings_config.append(w)
                    delimiter = \
                       misc.check_delimiter(sbtab_config.return_table_string())
                    (session.parameter_dict,
                     warnings) = misc.readout_config(sbtab_config)
                    if warnings != []:
                        for w in warnings:
                            session.warnings_config.append(w)
                except:
                    session.warnings_config.append('Error: The options SBtab f'
                                                   'ile could not be validated'
                                                   '. We suggest to validate t'
                                                   'he file manually on www.sb'
                                                   'tab.net.')
        except:
            session.warnings_config.append('Error: This is not a valid SBtab o'
                                           'ptions file.')

    elif configform.errors: response.flash = 'Form has errors'

    ###########################################################################
    ### RIGHT SIDE: Fast lane balancing (SBtab only)                       ####
    ###########################################################################

    sbtab_fl_form = SQLFORM.factory(Field('File', 'upload', uploadfolder="/tmp",
                                          label='Upload SBtab file (.tsv)',
                                          requires=IS_LENGTH(10485760, 1,
                                                             error_message='Max upload size: 10MB')))

    # file is uploaded, checked, and added to the list; variables are declared
    if sbtab_fl_form.process(formname='form_five').accepted:
        response.flash = 'form accepted'
        session.warnings_fl = []
        filename = request.vars.File.filename
        sbtab_file = request.vars.File.value
        session.warnings = []
        if 'sbtab_fls' not in session:
            session.sbtab_fls = []
            session.sbtab_fl_names = []

        if 'sbmls' not in session:
            session.sbmls = []
            session.sbml_names = []

        if 'sbtabs' not in session:
            session.sbtabs = []
            session.sbtab_names = []

        if 'priors' not in session:
            session.priors = []
            session.prior_names = []
            prior_open = open('./applications/pb/static/files/default_files/pb'
                              '_prior.tsv')
            prior_file = prior_open.read()
            sbtab_prior = SBtab.SBtabTable(prior_file, 'pb_prior.tsv')
            session.prior = sbtab_prior
            session.prior_name = 'pb_prior.tsv'
            session.priors.append(session.prior)
            session.prior_names.append(session.prior_name)

        try:
            sbtab_fl = SBtab.SBtabTable(sbtab_file.decode('utf-8'), filename)
            valid_extension = misc.validate_file_extension(sbtab_fl.filename,
                                                           'sbtab')
            if not valid_extension:
                session.warnings_fl.append('Error: The file %s has not a valid'
                                           ' file extension. Please upload a v'
                                           'alid SBtab file.' % (filename))
            if filename in session.sbtab_fl_names:
                session.warnings_fl.append('Error: Duplicate file name. Please'
                                           ' remove the file with the same nam'
                                           'e before uploading.')
            elif valid_extension:
                session.sbtab_fls.append(sbtab_file.decode('utf-8'))
                session.sbtab_fl_names.append(sbtab_fl.filename)
        except:
            session.warnings_fl.append('Error: The file %s is not compliant wi'
                                       'th SBtab. Please upload a valid SBtab '
                                       'file.' % (filename))

        try: redirect(URL('../default/balancing'))
        except: redirect(URL('../balancing'))
    elif sbtab_fl_form.errors:
        response.flash = 'Form has errors'

    ###########################################################################
    ###                       BUTTONS                                     #####
    ###########################################################################

    ### 1: SBML
    if request.vars.erase_button_sbml:
        del session.sbmls[int(request.vars.erase_button_sbml)]
        del session.sbml_names[int(request.vars.erase_button_sbml)]
        session.sbml = None
        session.sbml_name = None
        session.warnings_sbml = []
        try: redirect(URL('../default/balancing'))
        except: redirect(URL('../balancing'))

    ### 2. SBtab
    if request.vars.erase_button_sbtab:
        del session.sbtabs[int(request.vars.erase_button_sbtab)]
        del session.sbtab_names[int(request.vars.erase_button_sbtab)]
        session.sbtab = None
        session.sbtab_name = None
        session.warnings_sbtab = []
        try: redirect(URL('../default/balancing'))
        except: redirect(URL('../balancing'))

    ### 3: Prior
    if request.vars.erase_button_prior:
        del session.priors[int(request.vars.erase_button_prior)]
        del session.prior_names[int(request.vars.erase_button_prior)]
        session.prior = None
        session.prior_name = None
        session.warnings_prior = []
        try: redirect(URL('../default/balancing'))
        except: redirect(URL('../balancing'))

    ### 4: Config
    if request.vars.erase_button_config:
        session.config_file = None
        session.config_filename = None
        session.warnings_config = []
        session.parameter_dict = {}
        try: redirect(URL('../default/balancing'))
        except: redirect(URL('../balancing'))

    ### 5: Fast Lane
    if request.vars.erase_button_fl:
        del session.sbtab_fls[int(request.vars.erase_button_fl)]
        del session.sbtab_fl_names[int(request.vars.erase_button_fl)]
        session.sbtab_fl = None
        session.sbtab_fl_name = None
        session.warnings_fl = []
        try: redirect(URL('../default/balancing'))
        except: redirect(URL('../balancing'))

    ### 6: If stuff is in place, we can perform the classic balancing
    if request.vars.balance_classic:
        if 'priors' not in session:
            session.priors = []
            session.prior_names = []
            prior_open = open('./applications/pb/static/files/default_files/pb'
                              '_prior.tsv')
            prior_file = prior_open.read()
            sbtab_prior = SBtab.SBtabTable(prior_file, 'pb_prior.tsv')
            session.prior = sbtab_prior
            session.prior_name = 'pb_prior.tsv'
            session.priors.append(session.prior)
            session.prior_names.append(session.prior_name)

        # get SBML
        sbml_file = session.sbml
        sbml_filename = session.sbml_name
        # get SBtab or produce empty SBtab
        if 'sbtab' in session:
            if session.sbtab is not None:
                sbtab_file = session.sbtab
                sbtab_filename = session.sbtab_name
                session.emptysbtab = False
            else: session.emptysbtab = True
        else: session.emptysbtab = True

        # get prior file or open default
        if 'prior' in session:
            if session.prior is not None:
                sbtab_prior = session.prior
                prior_filename = session.prior_name
                emptyprior = False
            else: emptyprior = True
        else: emptyprior = True

        if emptyprior:
            try:
                prior_open = open('./applications/pb/static/files/default_file'
                                  's/pb_prior.tsv')
                prior_file = prior_open.read()
                sbtab_prior = SBtab.SBtabTable(prior_file, 'pb_prior.tsv')
                session.prior = sbtab_prior
                session.prior_name = sbtab_prior.filename
                session.priors.append(session.prior)
                session.prior_names.append(session.prior_name)
                emptyprior = False
            except:
                session.warnings_prior.append('Error loading the default prior'
                                              ' table.')

        # get config file if provided
        if 'config' in session:
            if session.config is not None:
                sbtab_config = session.config_file
                config_filename = session.config_filename
                if sbtab_config.table_type != 'PbConfig':
                    session.warnings_config.append('Error: The SBtab file has '
                                                   'an incorrect table type'
                                                   ': "'"%s"'". The correct ta'
                                                   'ble type would be "'"PbConfig"'".'
                                                   '\n' % (TableValidClass.sbtab.table_type))

                # get default definition file
                try:
                    def_file_open = open('./applications/pb/static/files/defau'
                                         'lt_files/definitions.tsv')
                    definition_file = def_file_open.read()
                    definition_name = 'definitions.tsv'
                    sbtab_def = SBtab.SBtabTable(definition_file, definition_name)
                    TableValidClass = validatorSBtab.ValidateTable(sbtab_config,
                                                                   sbtab_def)
                    warnings = TableValidClass.return_output()
                    if warnings != []:
                        for w in warnings:
                            session.warnings_config.append(itemx)
                    delimiter = \
                       misc.check_delimiter(sbtab_config.return_table_string())
                    (session.parameter_dict, warnings) = \
                                            misc.readout_config(sbtab_config)
                    if warnings != []:
                        for w in warnings:
                            session.warnings_config.append(itemx)
                except:
                    session.warnings_config.append('Error: The options SBtab f'
                                                   'ile could not be validated'
                                                   '. We suggest to validate t'
                                                   'he file manually on www.sb'
                                                   'tab.net.')
                    try: redirect(URL('../default/balancing'))
                    except: redirect(URL('../balancing'))

        if 'parameter_dict' not in session or \
           'config' not in session.parameter_dict:
            session.parameter_dict = {'config': False}

        # 1: extract priors and pseudos
        try:
            (pseudos,
             priors,
             pmin,
             pmax) = misc.extract_pseudos_priors(sbtab_prior)
        except:
            session.warnings_prior.append('Error: The default prior table %s c'
                                          'ould not be processed properly'
                                          '.' % (sbtab_prior.filename))
            try: redirect(URL('../default/balancing'))
            except: redirect(URL('../balancing'))

        # 2: makeSBtab (either from an empty SBtab or from a given SBtab)
        if session.emptysbtab:
            try:
                reader = libsbml.SBMLReader()
                sbml = reader.readSBMLFromString(sbml_file)
                sbml_model = sbml.getModel()
                pb = balancer.ParameterBalancing(sbml_model)
                sbtab_data = pb.make_empty_sbtab(pmin, pmax,
                                                 session.parameter_dict)
            except:
                session.warnings_sbml.append('Error: The SBML file %s could no'
                                             't be processed properly'
                                             '.' % (sbml_filename))
                try: redirect(URL('../default/balancing'))
                except: redirect(URL('../balancing'))
        else:
            try:
                reader = libsbml.SBMLReader()
                sbml = reader.readSBMLFromString(sbml_file)
                sbml_model = sbml.getModel()
                pb = balancer.ParameterBalancing(sbml_model)
                sbtab_data = pb.make_sbtab(sbtab_file, sbtab_file.filename,
                                           'All organisms', 43,
                                           pmin, pmax, session.parameter_dict)
            except:
                session.warnings_sbml.append('Error: The SBML file %s could no'
                                             't be processed properly'
                                             '.' % (sbml_filename))
                session.warnings_sbtab.append('Error: The SBtab file %s could no'
                                             't be processed properly'
                                             '.' % (sbtab_filename))
                try: redirect(URL('../default/balancing'))
                except: redirect(URL('../balancing'))

        # 2b: update prior information for model if alternate prior is used
        try:
            pb.get_parameter_information(sbtab_prior)
        except:
            session.warnings_prior.append('The alternate prior could not'
                                          ' be loaded.')
                       
        # 3: fill them in the SBtab file
        try:
            if 'use_pseudo_values' in session.parameter_dict and \
               session.parameter_dict['use_pseudo_values'] == 'True':
                sbtab_old = copy.deepcopy(sbtab_data)
                sbtab_new = pb.fill_sbtab(sbtab_old, pseudos, priors)
                pseudo_flag = 'pseudos'
            else:
                sbtab_new = pb.fill_sbtab(sbtab_data)
                pseudo_flag = 'no_pseudos'
        except:
                session.warnings_sbtab.append('Error: The SBtab data could not'
                                              'be processed correctly.')
                try: redirect(URL('../default/balancing'))
                except: redirect(URL('../balancing'))

        # 4: construct parameter dictionary
        if 'temperature' not in session.parameter_dict.keys():
            session.parameter_dict['temperature'] = 300
        if 'ph' not in session.parameter_dict.keys():
            session.parameter_dict['ph'] = 7

        types = ['standard chemical potential',
                 'catalytic rate constant geometric mean', 'Michaelis constant',
                 'activation constant', 'inhibitory constant', 'concentration',
                 'concentration of enzyme', 'equilibrium constant',
                 'substrate catalytic rate constant',
                 'forward maximal velocity', 'product catalytic rate constant',
                 'reverse maximal velocity', 'chemical potential',
                 'reaction affinity']

        for typ in types:
            session.parameter_dict[typ] = True

        # 5: BALANCE PARAMETERS
        try:
            (sbtab_final, mean_vector, mean_vector_inc, c_post,
             c_post_inc, r_matrix, shannon,
             log) = pb.make_balancing(sbtab_new, sbtab_data,
                                      pmin, pmax, session.parameter_dict)
            session.log = log
        except:
            session.warnings_sbml.append('Error: The balancing process was err'
                                         'oneous. Please check your input file'
                                         's for validity.')
            try: redirect(URL('../default/balancing'))
            except: redirect(URL('../balancing'))

        # 7: fill the model with the parameters and the kinetics
        try: param = session.parameter_dict['parametrisation']
        except: param = 'hal'
        try: prefac = session.parameter_dict['enzyme_prefactor']
        except: prefac = True
        try: inh = session.parameter_dict['default_inh']
        except: inh = 'complete_inh'
        try: act = session.parameter_dict['default_act']
        except: act = 'complete_act'
        try: overwrite = session.parameter_dict['overwrite_kinetics']
        except: overwrite = True

        try:
            kineticizer_cs = kineticizer.KineticizerCS(sbml_model, sbtab_final,
                                                       param, prefac, inh, act,
                                                       overwrite)
        except:
            session.warnings_sbtab.append('Error: The parameters and kinetics '
                                          'could not be written to the output '
                                          'model.')
            try: redirect(URL('../default/balancing'))
            except: redirect(URL('../balancing'))

        # 8: post production of SBtab and SBML formats
        # 8.1: SBML
        try:
            ### model_name
            model_name = str(sbml_filename)[:-4] + '_balanced_model.xml'
            sbml_code = '<?xml version="1.0" encoding="UTF-8"?>\n' + \
                        sbml_model.toSBML()
            session.result_sbml = [sbml_code]
            session.result_sbml_name = [model_name]
        except:
            session.warnings_sbml.append('Error: The new SBML model could not '
                                         'be produced.')
            try: redirect(URL('../default/balancing'))
            except: redirect(URL('../balancing'))

        # 8.2: SBtab
        try:
            sbtab_file_new = open(sbtab_final.filename + '.tsv', 'w')
            sbtab_file_new.write(sbtab_final.return_table_string())
            sbtab_file_new.close()
            session.result_sbtab = [sbtab_final]
            session.result_sbtab_name = [sbtab_final.filename[:-4] + \
                                         '_balanced_parameters.tsv']
        except:
            session.warnings_sbtab.append('Error: The new SBtab file could not'
                                          ' be produced.')
            try: redirect(URL('../default/balancing'))
            except: redirect(URL('../balancing'))

        # 8.3: SBtab (for all)
        # Fix this as soon as we have the SBtab Document Class
        # Currently, we only save a string with all the single SBtabs in it,
        # but we need a list of SBtab objects and a proper show-function for it
        try:
            document = sbml2sbtab.SBMLDocument(sbml_model, model_name)
            sbtab_string = ''
            (sbtab_all, warnings) = document.makeSBtabs()
            for sbtab in sbtab_all:
                sbtab_string += sbtab.return_table_string() + '\n\n'
            sbtab_string += sbtab_final.return_table_string() + '\n\n'
            if 'prior' in session:
                sbtab_string += session.prior.return_table_string() + '\n\n'
            if 'config_file' in session and session.config_file is not None:
                sbtab_string += session.config_file.return_table_string() + \
                                '\n\n'
            session.result_sbtab.append(sbtab_string)
            session.result_sbtab_name.append(sbml_filename[:-4] + \
                                             '_balanced_model.tsv')
        except:
            session.warnings_sbtab.append('Error: It was not possible to save '
                                          'the model and parameters in one SBt'
                                          'ab file for download.')
            try: redirect(URL('../default/balanced'))
            except: redirect(URL('../balanced'))

        # remove prior again to not have it shown initially on balancing screen
        prior_file = None
        prior_filename = None
        session.priors = []
        session.prior_names = []
        session.config_file = None
        session.config_filename = None
        try: redirect(URL('../default/balanced'))
        except: redirect(URL('../balanced'))

    if request.vars.clearsession:
        clearsession()

    ################################################### CLASSIC#

    # FASTLANE##################################################

    if request.vars.balance_fastlane:
        session.sbtab_fl = session.sbtab_fls[0]
        session.sbtab_fl_name = session.sbtab_fl_names[0]
        sbtab_data = False
        sbtab_q = False
        sbtab_c = False
        session.warnings = []

        # check out which sbtabs are available
        try:
            sbtabs = misc.cut_tabs(session.sbtab_fl)
            for sbtab in sbtabs:
                if 'Quantity' in sbtab and '!MathematicalType' in sbtab:
                    sbtab_prior = SBtab.SBtabTable(sbtab,
                                                   session.sbtab_fl_name[:-4] + \
                                                   '_%s.tsv' % 'prior')
                elif 'Quantity' in sbtab and '!Mean':
                    sbtab_data = SBtab.SBtabTable(sbtab,
                                                  session.sbtab_fl_name[:-4] + \
                                                  '_%s.tsv' % 'data')
                elif 'PbConfig' in sbtab:
                    sbtab_config = SBtab.SBtabTable(sbtab,
                                                    session.sbtab_fl_name[:-4] + \
                                                    '_%s.tsv' % 'config')
        except:
            session.warnings_fl.append('Error: Could not separate the bundled '
                                       'SBtab file. Please check syntax validi'
                                       'ty.')
            try: redirect(URL('../default/balancing'))
            except: redirect(URL('../balancing'))

        # convert sbtab to sbml
        try:
            Conversion_class = sbtab2sbml.SBtabDocument(session.sbtab_fl,
                                                        session.sbtab_fl_name)
            (sbml_file, warnings) = Conversion_class.makeSBML()
            session.sbmls.append(sbml_file)
            session.sbml_names.append(session.sbtab_fl_name[:-4] + '.xml')
            if sbml_file is False:
                print('The SBML file was not computed correctly from the bundl'
                      'ed SBtab')
            if warnings != []:
                for w in warnings:
                    session.warnings_fl.append(w)
        except:
            session.warnings_fl.append('Error: The SBtab file did not include '
                                       'a valid Reaction SBtab that could be c'
                                       'onverted to SBML.')
            try: redirect(URL('../default/balancing'))
            except: redirect(URL('../balancing'))

        # check out prior file if available
        try:
            validity = misc.valid_prior(sbtab_prior)
            if validity != []:
                for v in validity:
                    session.warnings_fl.append(v)
            session.priors.append(sbtab_prior)
            session.prior_names.append('embedded_prior.tsv')
            session.prior = sbtab_prior
            session.prior_name = 'embedded_prior.tsv'
        except:
            prior_open = open('./applications/pb/static/files/default_files/pb'
                              '_prior.tsv')
            prior_file = prior_open.read()
            session.prior_name = 'pb_prior.tsv'
            sbtab_prior = SBtab.SBtabTable(prior_file, session.prior_name)
            session.prior = sbtab_prior
            session.priors = []
            session.prior_names = []
            try:
                session.priors.append(sbtab_prior)
                session.prior_names.append(session.prior_name)
            except:
                session.warnings_fl.append('Error loading the default prior table.')
                try: redirect(URL('../default/balancing'))
                except: redirect(URL('../balancing'))

        # 0: config:
        if sbtab_config:
            session.config_file = sbtab_config
            session.config_filename = 'embedded_options.tsv'
            try:
                def_file_open = open('./applications/pb/static/files/default_f'
                                     'iles/definitions.tsv')
                definition_file = def_file_open.read()
                definition_name = 'definitions.tsv'
                sbtab_def = SBtab.SBtabTable(definition_file, definition_name)
                TableValidClass = validatorSBtab.ValidateTable(sbtab_config,
                                                               sbtab_def)
                if sbtab_config.table_type != 'PbConfig':
                    session.warnings_fl.append('Error: The SBtab file has an i'
                                               'ncorrect table type: "'"%s"'".'
                                               ' The correct table type would '
                                               'be "'"PbConfig"'".\n' % (sbtab_config.table_type))
                warnings = TableValidClass.return_output()
                if warnings != []:
                    for w in warnings:
                        session.warnings_fl.append(w)
                (session.parameter_dict,
                 warnings) = misc.readout_config(sbtab_config)
                if warnings != []:
                    for w in warnings:
                        session.warnings_fl.append(w)
            except:
                session.warnings_fl.append('Error: The options SBtab file coul'
                                           'd not be validated. We suggest to '
                                           'validate the file manually on www.'
                                           'sbtab.net.')
                try: redirect(URL('../default/balancing'))
                except: redirect(URL('../balancing'))

        if 'parameter_dict' not in session:
            session.parameter_dict = {'config': False}

        # 1: extract priors and pseudos
        try:
            (pseudos,
             priors, pmin,
             pmax) = misc.extract_pseudos_priors(session.prior)
        except:
            session.warnings_fl.append('Error: The prior table %s could not be'
                                       ' processed properly'
                                       '.' % (session.prior_name))
            try: redirect(URL('../default/balancing'))
            except: redirect(URL('../balancing'))

        # 2: makeSBtab (either from an empty SBtab or from a given SBtab)
        if not sbtab_data:
            try:
                reader = libsbml.SBMLReader()
                sbml = reader.readSBMLFromString(sbml_file)
                sbml_model = sbml.getModel()
                pb = balancer.ParameterBalancing(sbml_model)
                sbtab_data = pb.make_empty_sbtab(pmin, pmax,
                                                 session.parameter_dict)
            except:
                session.warnings_fl.append('Error: The model information of %s'
                                           ' could not be processed properly'
                                           '.' % (session.sbtab_fl_name))
                try:
                    redirect(URL('../default/balancing'))
                except:
                    redirect(URL('../balancing'))
        else:
            try:
                reader = libsbml.SBMLReader()
                sbml = reader.readSBMLFromString(sbml_file)
                sbml_model = sbml.getModel()
                pb = balancer.ParameterBalancing(sbml_model)
                sbtab_data = pb.make_sbtab(sbtab_data, sbtab_data.filename,
                                           'All organisms', 43,
                                           pmin, pmax, session.parameter_dict)
            except:
                session.warnings_fl.append('Error: The model information of %s'
                                           ' could not be processed properly'
                                           '.' % (session.sbtab_fl_name))
                try: redirect(URL('../default/balancing'))
                except: redirect(URL('../balancing'))

        # 2b: update prior information for model if alternate prior is used
        try:
            pb.get_parameter_information(sbtab_prior)
        except:
            session.warnings_prior.append('The alternate prior could not'
                                          ' be loaded.')
                                
        # 3: fill them in the SBtab file
        try: pseudo = bool(session.parameter_dict['use_pseudo_values'])
        except: pseudo = True

        try:
            if pseudo:
                sbtab_old = copy.deepcopy(sbtab_data)
                sbtab_new = pb.fill_sbtab(sbtab_old, pseudos, priors)
                pseudo_flag = 'pseudos'
            else:
                sbtab_new = pb.fill_sbtab(sbtab_data)
                pseudo_flag = 'no_pseudos'
        except:
            session.warnings_fl.append('Error: The prior table %s could not be'
                                       ' processed properly'
                                       '.' % (session.prior_name))
            try: redirect(URL('../default/balancing'))
            except: redirect(URL('../balancing'))

        # 4: construct parameter dictionary
        if 'temperature' not in session.parameter_dict.keys():
            session.parameter_dict['temperature'] = 300
        else: session.parameter_dict['temperature'] = \
            float(session.parameter_dict['temperature'])
        if 'ph' not in session.parameter_dict.keys():
            session.parameter_dict['ph'] = 7
        else:
            session.parameter_dict['ph'] = float(session.parameter_dict['ph'])

        types = ['standard chemical potential',
                 'catalytic rate constant geometric mean',
                 'Michaelis constant', 'activation constant',
                 'inhibitory constant', 'concentration',
                 'concentration of enzyme', 'equilibrium constant',
                 'substrate catalytic rate constant',
                 'forward maximal velocity', 'product catalytic rate constant',
                 'reverse maximal velocity', 'chemical potential',
                 'reaction affinity']

        for typ in types: session.parameter_dict[typ] = True

        #####
        # 5: BALANCE PARAMETERS
        try:
            (sbtab_final, mean_vector, mean_vector_inc, c_post, c_post_inc,
             r_matrix, shannon, log) = pb.make_balancing(sbtab_new,
                                                         sbtab_data,
                                                         pmin, pmax,
                                                         session.parameter_dict)
            session.log = log
        except:
            session.warnings_fl.append('Error: The balancing process was erron'
                                       'eous. Please check validity of input f'
                                       'iles.')
            try: redirect(URL('../default/balancing'))
            except: redirect(URL('../balancing'))

        # 7: fill the model with the parameters and the kinetics
        try:
            try: param = session.parameter_dict['parametrisation']
            except: param = 'hal'
            try: prefac = session.parameter_dict['enzyme_prefactor']
            except: prefac = True
            try: inh = session.parameter_dict['default_inh']
            except: inh = 'complete_inh'
            try: act = session.parameter_dict['default_act']
            except: act = 'complete_act'
            try: overwrite = session.parameter_dict['overwrite_kinetics']
            except: overwrite = True
            kineticizer_cs = kineticizer.KineticizerCS(sbml_model, sbtab_final,
                                                       param, prefac, inh, act,
                                                       overwrite)
        except:
            session.warnings_fl.append('Error: The parameters and kinetics cou'
                                       'ld not be written to the model.')
            try: redirect(URL('../default/balancing'))
            except: redirect(URL('../balancing'))

        # 8.1: Output SBML
        try:
            ### model_name
            model_name = str(session.sbtab_fl_name)[:-4] + \
                         '_balanced_model.xml'
            sbml_code = '<?xml version="1.0" encoding="UTF-8"?>\n' + \
                        sbml_model.toSBML()
            session.result_sbml = [sbml_code]
            session.result_sbml_name = [model_name]
        except:
            session.warnings_fl.append('Error: The new SBML model could not be'
                                       ' produced.')
            try: redirect(URL('../default/balancing'))
            except: redirect(URL('../balancing'))

        # 8.2: Output SBtab
        try:
            session.result_sbtab = [sbtab_final]
            session.result_sbtab_name = [session.sbtab_fl_name[:-4] + \
                                         '_balanced_parameters.tsv']
        except:
            session.warnings_fl.append('Error: The new SBtab file could not be'
                                       ' produced.')
            try: redirect(URL('../default/balancing'))
            except: redirect(URL('../balancing'))

        # 8.3: Output SBtab (for all)
        try:
            sbtabs_all = ''
            for sbtab in sbtabs:
                sbtabs_all += sbtab + '\n'
            sbtabs_all += sbtab_final.return_table_string() + '\n'
            session.result_sbtab.append(sbtabs_all)
            session.result_sbtab_name.append(session.sbtab_fl_name[:-4] + \
                                             '_balanced_model.tsv')
            session.result_sbtab_string = sbtabs_all
            session.result_sbtab_string_name = session.sbtab_fl_name[:-4] + \
                                               '_balanced_model.tsv'
        except:
            session.warnings_fl.append('Error: It was not possible to save the'
                                       ' model and parameters in one SBtab fil'
                                       'e for download.')
            try: redirect(URL('../default/balancing'))
            except: redirect(URL('../balancing'))

        # remove prior and config again to not have it shown initially on
        # balancing screen
        prior_file = None
        prior_filename = None
        session.priors = []
        session.prior_names = []
        session.config_file = None
        session.config_filename = None
        try: redirect(URL('../default/balanced'))
        except: redirect(URL('../balanced'))

        ################################## BALANCING

    return dict(SBML_FORM=sbmlform, SBMLS=session.sbmls,
                SBML_NAMES=session.sbml_names,
                WARNINGS_SBML=session.warnings_sbml, SBML=session.sbml,
                SBML_NAME=session.sbml_name, SBTAB_FORM=sbtabform,
                SBTABS=session.sbtabs, SBTAB_NAMES=session.sbtab_names,
                WARNINGS_SBTAB=session.warnings_sbtab,
                WARNINGS_SBTAB_FL=session.warnings_fl, SBTAB=session.sbtab,
                SBTAB_NAME=session.sbtab_name, PRIOR_FORM=priorform,
                WARNINGS_PRIOR=session.warnings_prior, PRIORS=session.priors,
                PRIOR_NAMES=session.prior_names, PRIOR=session.prior,
                PRIOR_NAME=session.prior_name,
                WARNINGS_CONFIG=session.warnings_config,
                CONFIG_FILE=session.config_file,
                CONFIG_FILENAME=session.config_filename,
                CONFIG_FORM=configform, SBTAB_FL_FORM=sbtab_fl_form,
                SBTAB_FLS=session.sbtab_fls,
                SBTAB_FL_NAMES=session.sbtab_fl_names,
                WARNINGS=session.warnings,
                RESULT_SBML_NAME=session.result_sbml_name,
                RESULT_SBML=session.result_sbml,
                RESULT_SBTAB_NAME=session.result_sbtab_name,
                RESULT_SBTAB=session.result_sbtab,
                RESULT_SBTAB_STRING=session.result_sbtab_string,
                RESULT_SBTAB_STRING_NAME=session.result_sbtab_string_name,
                LOG_FILE=session.log_file)


def balanced():
    '''
    fourth and last page of online balancing:
    view and/or download the output files
    '''
    response.title = T('Parameter Balancing for Kinetic Models of Cell '
                       'Metabolism')
    response.subtitle = T('Online Balancing')

    if request.vars.download_button_sbml:
        download_sbml()
        pass

    if request.vars.erase_button_sbml2:
        del session.result_sbml[int(request.vars.erase_button_sbml2)]
        del session.result_sbml_name[int(request.vars.erase_button_sbml2)]
        session.log = False
        try: redirect(URL('../default/balanced'))
        except: redirect(URL('../balanced'))

    if request.vars.download_button_sbtab:
        download_sbtab()
        pass

    if request.vars.download_log_button:
        download_log()
        pass

    if request.vars.erase_button_sbtab2:
        del session.result_sbtab[int(request.vars.erase_button_sbtab2)]
        del session.result_sbtab_name[int(request.vars.erase_button_sbtab2)]
        session.warnings_sbtab = []
        try: redirect(URL('../default/balanced'))
        except: redirect(URL('../balanced'))

    return dict(SBML_NAME=session.sbml_name, SBMLS=session.sbmls,
                SBML_NAMES=session.sbml_names, SBTABS=session.sbtabs,
                SBTAB_NAMES=session.sbtab_names, LOG_FILE=session.log,
                WARNINGS=session.warnings, RESULT_SBML=session.result_sbml,
                RESULT_SBTAB_NAME=session.result_sbtab_name,
                RESULT_SBTAB=session.result_sbtab,
                RESULT_SBTAB_STRING_NAME=session.result_sbtab_string_name,
                RESULT_SBTAB_STRING=session.result_sbtab_string,
                RESULT_SBML_NAME=session.result_sbml_name)


def show_sbml():
    '''
    function that converts SBML file to HTML in order to display it in browser
    '''
    return misc.xml2html(str(session.sbmls[int(request.args(0))]))


def show_sbml2():
    '''
    function that converts SBML file to HTML in order to display it in browser
    '''
    return misc.xml2html(session.result_sbml[int(request.args(0))])


def show_sbtab():
    '''
    function that converts SBtab file to HTML in order to display it in browser
    '''
    try:
        sbtab = session.sbtabs[int(request.args(0))]
    except: return 'The requested SBtab file cannot be loaded.'

    try:
        return misc.tsv_to_html(sbtab)
    except: return 'The requested SBtab file cannot be displayed.'


def show_sbtab2():
    '''
    function that converts SBtab file to HTML in order to display it in browser
    '''
    try:
        sbtab = session.result_sbtab[int(request.args(0))]
        filename = session.result_sbtab_name[int(request.args(0))]
    except: return 'The requested SBtab file cannot be loaded.'

    try:
        return misc.tsv_to_html(sbtab, filename)
    except: return 'The requested SBtab file cannot be displayed.'


def show_log():
    '''
    function that converts log file to HTML in order to display it in browser
    '''
    try:
        log_name = 'log_file.txt'
        log_file = session.log
    except: return 'The requested log file cannot be loaded.'

    try:
        delimiter = '\t'
        return misc.tsv_to_html(log_file, log_name)
    except: return 'The requested log file cannot be displayed.'


def show_sbtab_fl():
    '''
    function that converts fast lane SBtab file to HTML
    in order to display it in the browser
    '''
    try:
        file_name = session.sbtab_fl_names[int(request.args(0))]
    except: return 'The requested SBtab file cannot be loaded.'

    # This needs to be updated when the SBtab Document class is at hand
    try:
        return misc.tsv_to_html(sbtab_file, file_name)
    except: return 'The requested SBtab file cannot be displayed.'


def show_prior():
    '''
    function that converts prior SBtab file to HTML
    in order to display it in the browser
    '''
    try:
        sbtab_prior = session.priors[int(request.args(0))]
    except: return 'The requested prior table cannot be loaded.'

    try:
        return misc.tsv_to_html(sbtab_prior)
    except: return 'The requested prior file cannot be displayed.'


def show_sbtab_conf():
    '''
    function that converts configure SBtab file to HTML
    in order to display it in the browser
    '''
    try:
        sbtab_config = session.config_file
    except: return 'The requested config table cannot be loaded.'

    try:
        return misc.tsv_to_html(sbtab_config)
    except: return 'The requested options file cannot be displayed.'


def download_sbtab():
    '''
    function for download of SBtab files
    '''
    response.headers['Content-Type'] = 'text/csv'
    attachment = 'attachment;filename=' + \
                 session.result_sbtab_name[int(request.vars.download_button_sbtab)]
    response.headers['Content-Disposition'] = attachment

    # fix the except part: we need something new, when we have
    # the SBtab Document class
    try:
        content = \
          session.result_sbtab[int(request.vars.download_button_sbtab)].return_table_string()
    except:
        content = \
          session.result_sbtab[int(request.vars.download_button_sbtab)]

    raise HTTP(200, str(content),
               **{'Content-Type': 'text/csv',
                  'Content-Disposition': attachment + ';'})


def download_sbml():
    '''
    function for download of SBML files
    '''
    response.headers['Content-Type'] = 'text/xml'
    attachment = 'attachment;filename=' + \
                 session.sbml_names[int(request.vars.download_button_sbml)]
    response.headers['Content-Disposition'] = attachment
    content = session.sbmls[int(request.vars.download_button_sbml)]

    raise HTTP(200, str(content),
               **{'Content-Type': 'text/xml',
                  'Content-Disposition': attachment + ';'})


def download_log():
    '''
    function for download of log files
    '''
    response.headers['Content-Type'] = 'text/csv'
    attachment = 'attachment;filename=' + session.sbml_name[:-4] + \
                 '_balancing_log.txt'
    response.headers['Content-Disposition'] = attachment
    content = session.log

    raise HTTP(200, str(content),
               **{'Content-Type': 'text/csv',
                  'Content-Disposition': attachment + ';'})


def user():
    """
    exposes:
    http://..../[app]/default/user/login
    http://..../[app]/default/user/logout
    http://..../[app]/default/user/register
    http://..../[app]/default/user/profile
    http://..../[app]/default/user/retrieve_password
    http://..../[app]/default/user/change_password
    http://..../[app]/default/user/manage_users (requires membership in
    use @auth.requires_login()
        @auth.requires_membership('group name')
        @auth.requires_permission('read', 'table name', record_id)
    to decorate functions that need access control
    """
    return dict(form=auth())


@cache.action()
def download():
    """
    allows downloading of uploaded files
    http://..../[app]/default/download/[filename]
    """
    return response.download(request, db)


def call():
    """
    exposes services. for example:
    http://..../[app]/default/call/jsonrpc
    decorate with @services.jsonrpc the functions to expose
    supports xml, json, xmlrpc, jsonrpc, amfrpc, rss, csv
    """
    return service()


@auth.requires_signature()
def data():
    """
    http://..../[app]/default/data/tables
    http://..../[app]/default/data/create/[table]
    http://..../[app]/default/data/read/[table]/[id]
    http://..../[app]/default/data/update/[table]/[id]
    http://..../[app]/default/data/delete/[table]/[id]
    http://..../[app]/default/data/select/[table]
    http://..../[app]/default/data/search/[table]
    but URLs must be signed, i.e. linked with
      A('table',_href=URL('data/tables',user_signature=True))
    or with the signed load operator
      LOAD('default','data.load',args='tables',ajax=True,user_signature=True)
    """
    return dict(form=crud())
