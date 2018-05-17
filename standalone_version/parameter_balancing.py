#!/usr/bin/env python
import copy
import libsbml
import os
import re
import sys
import argparse
try:
    from . import balancer
    from . import kineticizer
    from . import misc
    from . import SBtab
    from . import validatorSBtab
except:
    import balancer
    import kineticizer
    import misc
    import SBtab
    import validatorSBtab


def parameter_balancing_wrapper(parser_args):
    '''
    wrapper for the parameter balancing via the command line.
    receives the name of the SBML model and optional file names
    of SBtabs that need to be determined
    '''
    model_name = args.sbml
    parameter_dict = {}
    log_file = 'Parameter balancing log file of model %s\n' % (model_name)
    warn_flag = False

    ###########################
    # 1: open and prepare the files; then check for some rudimentary validity:
    # 1.1: SBML model
    reader = libsbml.SBMLReader()
    try: sbml_file = reader.readSBML(model_name)
    except:
        print('The SBML file %s could not be found.' % (model_name))
        sys.exit()
    try: sbml_model = sbml_file.getModel()
    except:
        print('The SBML file %s is corrupt. I quit.' % (model_name))
        sys.exit()
    valid_extension = misc.validate_file_extension(model_name, 'sbml')
    if not valid_extension:
        print('The SBML file %s has not the correct xml extension'\
              '. I quit.''' % (model_name))
        sys.exit()

    pb = balancer.ParameterBalancing(sbml_model)

    ###########################
    # 1.2: open and prepare the optional SBtab data file
    if args.sbtab_data:
        valid_extension = misc.validate_file_extension(args.sbtab_data,
                                                       'sbtab')
        if not valid_extension:
            print('The SBtab data file %s has not the correct file '\
                  'extension.' % (args.sbtab_data))
        try:
            f = open(args.sbtab_data, 'r')
            f_content = f.read()
        except:
            print('The SBtab data file %s cannot be found or'\
                  'read.' % args.sbtab_data)
        try: sbtab_delimiter = misc.check_delimiter(f_content)
        except: sbtab_delimiter = '\t'

        sbtab_data = SBtab.SBtabTable(f_content, args.sbtab_data)
        sbtab_data_validate = validatorSBtab.ValidateTable(sbtab_data,
                                                           args.sbtab_data)
        warnings = sbtab_data_validate.return_output()
        if warnings != []:
            warn_flag = True
            log_file += 'Log warnings for SBtab data file: '\
                        '%s\n' % args.sbtab_data
            for warning in warnings:
                log_file += warning + '\n'

    ###########################
    # 1.3: open and prepare an optional SBtab prior file;
    #      if this is not provided, open the default prior file
    if args.sbtab_prior:
        # try to open and read the file
        valid_extension = misc.validate_file_extension(args.sbtab_prior,
                                                       'sbtab')
        if not valid_extension:
            print('The SBtab prior file %s has not the correct file'\
                  'extension.' % (args.sbtab_prior))
        try:
            f = open(args.sbtab_prior, 'r')
            f_content = f.read()
        except:
            print('The SBtab prior file %s cannot be found or'\
                  'read.' % args.sbtab_prior)

        # initialise an SBtab object with the content and check its validity
        sbtab_prior = SBtab.SBtabTable(f_content, args.sbtab_prior)
        sbtab_prior_validate = validatorSBtab.ValidateTable(sbtab_prior,
                                                            args.sbtab_prior)

        # register warnings
        warnings = sbtab_prior_validate.return_output()
        if warnings != []:
            warn_flag = True
            log_file += 'Log warnings for SBtab prior file: '\
                        '%s\n\n' % args.sbtab_prior
            for warning in warnings:
                log_file += warning + '\n'

        valid_prior = misc.valid_prior(sbtab_prior)
        if valid_prior != []:
            warn_flag = True
            log_file += 'Log warnings for SBtab prior file: '\
                        '%s\n\n' % args.sbtab_prior
            for element in valid_prior:
                log_file += str(element) + '\n'

        # extract crucial information from prior
        (pseudos, priors, pmin, pmax) = misc.extract_pseudos_priors(sbtab_prior)
    else:
        # open default prior file
        p = os.path.dirname(os.path.abspath(__file__)) + '/files/default_'\
            'files/pb_prior.tsv'
        try:
            prior_file = open(p, 'r')
            prior = prior_file.read()
        except:
            print('The prior file (/files/default_files/pb_prior.tsv) coul'\
                  'd not be found. I quit.')
            sys.exit()

        sbtab_prior = SBtab.SBtabTable(prior, 'pb_prior.tsv')
        sbtab_prior_validate = validatorSBtab.ValidateTable(sbtab_prior,
                                                            'pb_prior.tsv')

        # register warnings
        warnings = sbtab_prior_validate.return_output()
        if warnings != []:
            warn_flag = True
            log_file += 'Log warnings for SBtab prior file: '\
                        '%s\n\n' % args.sbtab_prior
            for warning in warnings:
                log_file += warning + '\n'

        valid_prior = misc.valid_prior(sbtab_prior)
        if valid_prior != []:
            warn_flag = True
            log_file += 'Log warnings for SBtab prior file: '\
                        '%s\n\n' % args.sbtab_prior
            for element in valid_prior:
                log_file += str(element) + '\n'

        # extract crucial information from prior
        (pseudos, priors, pmin, pmax) = misc.extract_pseudos_priors(sbtab_prior)

    ###########################
    # 1.4: open and prepare an optional SBtab options file;
    #      if this is not provided, open the default options file
    if args.sbtab_options:
        valid_extension = misc.validate_file_extension(args.sbtab_options, 'sbtab')
        if not valid_extension:
            print('The SBtab options file %s has not the correct file'\
                  ' extension.' % (args.sbtab_options))
        try:
            f = open(args.sbtab_options, 'r')
            f_content = f.read()
        except:
            print('The SBtab options file %s cannot be found or'\
                  'read.' % args.sbtab_options)

        sbtab_options = SBtab.SBtabTable(f_content, args.sbtab_options)
        sbtab_options_validate = validatorSBtab.ValidateTable(sbtab_options,
                                                              args.sbtab_options)

        # register warnings
        warnings = sbtab_options_validate.return_output()
        if warnings != []:
            warn_flag = True
            log_file += 'Log warnings for SBtab options file: '\
                        '%s\n\n' % args.sbtab_options
            for warning in warnings:
                log_file += warning + '\n'

        (parameter_dict, log) = misc.readout_config(sbtab_options)
        if log != []:
            warn_flag = True
            log_file += 'Log warnings for SBtab options file: '\
                        '%s\n\n' % args.sbtab_options
            for element in log:
                log_file += str(element) + '\n'
    else:
        o = os.path.dirname(os.path.abspath(__file__)) + '/files/default_'\
            'files/pb_options.tsv'
        try:
            options_file = open(o, 'r')
            f_content = options_file.read()
        except:
            print('The options file (/files/default_files/pb_options.tsv) coul'\
                  'd not be found. I quit.')
            sys.exit()

        sbtab_options = SBtab.SBtabTable(f_content, 'pb_options.tsv')
        sbtab_options_validate = validatorSBtab.ValidateTable(sbtab_options,
                                                              'pb_options.tsv')

        # register warnings
        warnings = sbtab_options_validate.return_output()
        if warnings != []:
            warn_flag = True
            log_file += 'Log warnings for SBtab options file: '\
                        'pb_options.tsv\n\n'
            for warning in warnings:
                log_file += warning + '\n'

        (parameter_dict, log) = misc.readout_config(sbtab_options)

        if log != []:
            warn_flag = True
            log_file += 'Log warnings for SBtab options file: '\
                        '%s\n\n' % args.sbtab_options
            for element in log:
                log_file += str(element) + '\n'

    # Make empty SBtab if required
    if args.sbtab_data:
        sbtab = pb.make_sbtab(sbtab_data, args.sbtab_data, 'All organisms', 43,
                              pmin, pmax, parameter_dict)
        sbtabid2sbmlid = misc.id_checker(sbtab, sbml_model)
        if sbtabid2sbmlid != []:
            warn_flag = True
            log_file += 'Log warnings for SBtab data file: '\
                        '%s\n\n' % args.sbtab_data
            for element in sbtabid2sbmlid:
                log_file += element + '\n'
    else:
        sbtab = pb.make_empty_sbtab(pmin, pmax, parameter_dict)

    # end of file read in and processing;
    # now verify that all required keys are given for the parameter_dict;
    # if not, add them
    if 'temperature' not in parameter_dict.keys():
        parameter_dict['temperature'] = 300
    if 'ph' not in parameter_dict.keys(): parameter_dict['ph'] = 7
    if 'standard chemical potential' not in parameter_dict.keys():
        parameter_dict['standard chemical potential'] = True
    if 'catalytic rate constant geometric mean' not in parameter_dict.keys():
        parameter_dict['catalytic rate constant geometric mean'] = True
    if 'Michaelis constant' not in parameter_dict.keys():
        parameter_dict['Michaelis constant'] = True
    if 'activation constant' not in parameter_dict.keys():
        parameter_dict['activation constant'] = True
    if 'inhibitory constant' not in parameter_dict.keys():
        parameter_dict['inhibitory constant'] = True
    if 'concentration' not in parameter_dict.keys():
        parameter_dict['concentration'] = True
    if 'concentration of enzyme' not in parameter_dict.keys():
        parameter_dict['concentration of enzyme'] = True
    if 'equilibrium constant' not in parameter_dict.keys():
        parameter_dict['equilibrium constant'] = True
    if 'substrate catalytic rate constant' not in parameter_dict.keys():
        parameter_dict['substrate catalytic rate constant'] = True
    if 'product catalytic rate constant' not in parameter_dict.keys():
        parameter_dict['product catalytic rate constant'] = True
    if 'forward maximal velocity' not in parameter_dict.keys():
        parameter_dict['forward maximal velocity'] = True
    if 'reverse maximal velocity' not in parameter_dict.keys():
        parameter_dict['reverse maximal velocity'] = True
    if 'chemical potential' not in parameter_dict.keys():
        parameter_dict['chemical potential'] = True
    if 'reaction affinity' not in parameter_dict.keys():
        parameter_dict['reaction affinity'] = True
    if 'use_pseudo_values' not in parameter_dict.keys():
        parameter_dict['use_pseudo_values'] = False

    if args.verbose:
        print('\nFiles successfully read. Start balancing.')

    # 2: Parameter balancing
    if parameter_dict['use_pseudo_values'] == 'True' and not args.pb_pseudos:
        sbtab_old = copy.deepcopy(sbtab)
        sbtab_new = pb.fill_sbtab(sbtab_old, pseudos, priors)
        pseudo_flag = 'pseudos'
        if args.verbose:
            print('Parameter balancing is using pseudo values.')
    else:
        sbtab_new = pb.fill_sbtab(sbtab)
        pseudo_flag = 'no_pseudos'
        if args.verbose:
            print('Parameter balancing is not using pseudo values.')

    (sbtab_final, mean_vector, mean_vector_inc, c_post, c_post_inc,
     r_matrix, shannon, log) = pb.make_balancing(sbtab_new,
                                                 sbtab, pmin,
                                                 pmax,
                                                 parameter_dict)

    # 3: inserting parameters and kinetics into SBML model
    transfer_mode = {'standard chemical potential': 'weg',
                     'equilibrium constant': 'hal',
                     'catalytic rate constant': 'cat'}

    try:
        clear_mode = parameter_dict['parametrisation']
        mode = transfer_mode[clear_mode]
    except: mode = 'hal'
    try: enzyme_prefac = parameter_dict['prefac']
    except: enzyme_prefac = True
    try: def_inh = parameter_dict['default_inh']
    except: def_inh = 'complete_inh'
    try: def_act = parameter_dict['default_act']
    except: def_act = 'complete_act'
    try: overwrite = parameter_dict['overwrite']
    except: overwrite = True
    kineticizer_cs = kineticizer.KineticizerCS(sbml_model, sbtab_final, mode,
                                               enzyme_prefac, def_inh,
                                               def_act, True)

    if args.output_name:
        output_name = args.output_name
    else:
        try: rm = re.match('.*/(.*)', str(model_name)).group(1)[:-4]
        except: rm = str(model_name)[:-4]
        output_name = rm + '_balanced'

    if args.verbose:
        print('Done... writing output files.')

    if warn_flag:
        print('The parameter balancing issued warnings. Please generate the '
        'log file with the -l flag and check the warnings.')
        
    # 5: If requested write log file
    if args.pb_log:
        if log_file.count('\n') == 2:
            log_file += 'No warnings detected. \n'
        log = open(output_name + '_log.txt', 'w')
        log.write(log_file)
        log.close()
        if args.verbose:
            print('The log file %s has been written.' % (output_name + '_log.txt'))

    # 6: Write SBtab and SBML model
    sbtab_file_new = open(output_name + '.tsv', 'w')
    sbtab_file_new.write(sbtab_final.return_table_string())
    sbtab_file_new.close()
    if args.verbose:
        print('The SBtab file %s has been written.' % (output_name + '.tsv'))

    sbml_code = '<?xml version="1.0" encoding="UTF-8"?>\n' + sbml_model.toSBML()
    sbml_model_new = open(output_name + '.xml', 'w')
    sbml_model_new.write(sbml_code)
    sbml_model_new.close()
    if args.verbose:
        print('The SBML file %s has been written.' % (output_name + '.xml'))
        print('>> Goodbye.')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('sbml', help='Path to an SBML file.')
    parser.add_argument('--sbtab_data', help='Path to an SBtab data file.')
    parser.add_argument('--sbtab_prior', help='Path to an SBtab prior file.')
    parser.add_argument('--sbtab_options', help='Path to an SBtab options file.')
    parser.add_argument('--output_name', help='Choose a name for the output files.')
    parser.add_argument('-l', '--pb_log', help='Flag to print a log file.', action='store_true')
    parser.add_argument('-p', '--pb_pseudos', help='Flag for usage of pseudo values.', action='store_true')
    parser.add_argument('-v', '--verbose', help='Flag to display script messages.', action='store_true')

    args = parser.parse_args()

    parameter_balancing_wrapper(args)

