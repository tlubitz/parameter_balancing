#!/usr/bin/env python
import libsbml
import os
import re
import sys
from . import balancer
from . import kineticizer
from . import misc
from . import SBtab

def parameter_balancing_wrapper(model_name, first = None,
                                second = None, third = None):
    '''
    wrapper for the parameter balancing via the command line.
    receives the name of the SBML model and optional file names
    of SBtabs that need to be determined
    '''
    log_file = 'Parameter balancing log file of model %s\n\n' % (model_name)

    # 1: open and prepare the files; then check for some rudimentary validity:
    # 1.1: SBML model
    reader = libsbml.SBMLReader()
    try: sbml_file = reader.readSBML(model_name)
    except: print('The SBML file %s could not be found.' % (sbml_name))
    try: sbml_model = sbml_file.getModel()
    except:
        print('The SBML file %s is corrupt. I quit.' % (sbml_name))
        sys.exit()
    valid_extension = misc.validate_file_extension(model_name, 'sbml')
    if not valid_extension:
        print('The SBML file %s has not the correct '\
              'xml extension. I quit.' % (model_name))
        sys.exit()

    pb = balancer.ParameterBalancing(sbml_model)

    # 1.2: If more arguments than the SBML file name are given, determine
    # what SBtab files we have here
    def assign_sbtab(sb, table_type, s):
        '''
        assigns the sbtab to the correct table type
        '''
        global sbtab
        global sbtab_name
        global prior
        global config
        if table_type == 'Quantity':
            sbtab = sb
            sbtab_name = s
        elif table_type == 'QuantityInfo': prior = sb
        elif table_type == 'PbConfig': config = sb
        else: print('The provided SBtab file %s could not be assigned '\
                    'properly' % (table_type))

    global sbtab
    global sbtab_name
    global prior
    global config
    sbtab = False
    sbtab_name = False
    prior = False
    config = False
    sbtabs = [first, second, third]

    for s in sbtabs:
        if s is not None:
            try:
                valid_extension = misc.validate_file_extension(s, 'sbtab')
                if not valid_extension:
                    print('The SBtab file %s has not the correct tsv '\
                          'extension.' % (first))
                undetermined = open(s, 'r')
                sb = undetermined.read()
                table_type = misc.table_type(sb)
                assign_sbtab(sb, table_type, s)
            except:
                print('The file %s could not be read and can thus not be'\
                      'integrated.' % (first))

    if sbtab_name is False: sbtab_name = 'output.csv'

    # PARAMETER FILE
    if sbtab:
        try: sbtab_delimiter = misc.check_delimiter(sb)
        except: sbtab_delimiter = '\t'
        sbtab = SBtab.SBtabTable(sbtab)
        no_sbtab = False
    else: no_sbtab = True

    # PRIOR FILE
    if prior:
        try: prior_delimiter = misc.check_delimiter(prior)
        except: prior_delimiter = '\t'
        valid_prior = misc.valid_prior(prior, prior_delimiter)
        if valid_prior != []:
            for element in valid_prior:
                log_file += str(element) + '\n'
    else:
        p = os.path.dirname(os.path.abspath(__file__)) + '/files/default_'\
            'files/pb_prior.tsv'
        try: prior_file = open(p, 'r')
        except:
            print('The prior file (/files/default_files/pb_prior.tsv) coul'\
                  'd not be found. I quit.')
            sys.exit()
        prior = prior_file.read()
        prior_delimiter = '\t'
    (pseudos, pmin, pmax) = misc.extract_pseudos(prior,
                                                 prior_delimiter)

    # CONFIG FILE
    if config:
        config_delimiter = misc.check_delimiter(config)
        (parameter_dict, log) = misc.readout_config(config, config_delimiter)
        if log != []:
            for element in log:
                log_file += str(element) + '\n'
    else:
        try:
            c = os.path.dirname(os.path.abspath(__file__)) + '/files/default_'\
                'files/pb_options.tsv'
            cf = open(c, 'r')
            config = cf.read()
            config_delimiter = '\t'
            (parameter_dict,log) = misc.readout_config(config, config_delimiter)
        except:
            print('The config file (/files/default_files/pb_config.tsv) coul'\
                  'd not be found and not config file was provided by user; '\
                  'the values are set to default.')
            parameter_dict = {}
    
    # MAKE EMPTY SBTAB IF REQUIRED
    if no_sbtab:
        sbtab = pb.make_empty_sbtab(pmin, pmax, parameter_dict)
    else:
        sbtab = pb.make_empty_sbtab(sbtab, sbtab_name, 'All organisms', 43, pmin, pmax, parameter_dict)
        sbtabid2sbmlid = misc.id_checker(sbtab, sbml_model)
        if sbtabid2sbmlid != []:
            for element in sbtabid2sbmlid:
                log_file.append(element)
               
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
       
    print('\nFiles successfully read. Start balancing.\n')

    # 2: Parameter balancing
    if parameter_dict['use_pseudo_values']:
        sbtab_new = pb.fill_sbtab(sbtab, pseudos)
        pseudo_flag = 'pseudos'
    else:
        sbtab_new = pb.fill_sbtab(sbtab)
        pseudo_flag = 'no_pseudos'
        
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
    kineticizer_cs = kineticizer.kineticizer_cs(sbml_model, sbtab_final, mode,
                                                enzyme_prefac, def_inh,
                                                def_act, True)

    output_name = input('\nEnter optional name for output files. If you do '\
                        'not enter a name, the files are named after the '\
                        'input file.\n')
    if output_name == '':
        try: rm = re.match('.*/(.*)', str(model_name)).group(1)[:-4]
        except: rm = str(model_name)[:-4]
        output_name = rm + '_balanced'

    print('\nDone... now writing output files to current working directory. '\
          '\n\nGoodbye!\n\n')

    # 5: Write SBtab and SBML model
    model_name = str(model_name)[:-4]
    sbtab_file_new = open(output_name + '.tsv', 'w')
    sbtab_file_new.write(sbtab_final.return_table_string())
    sbtab_file_new.close()
    
    sbml_code = '<?xml version="1.0" encoding="UTF-8"?>\n' + sbml_model.toSBML()
    sbml_model_new = open(output_name + '.xml', 'w')
    sbml_model_new.write(sbml_code)
    sbml_model_new.close()

    try:
        os.remove('cpost.txt')
        os.remove('medians.txt')
    except: pass

if __name__ == '__main__':

    try: model_name = sys.argv[1]
    except:
        print('\nPlease provide at least an SBML model as argument for the '\
              'script, i.e. \n >python parameter_balancing.py '\
              'your_sbml_file.xml \n')
        sys.exit()

    try:
        first = sys.argv[2]
        try:
            second = sys.argv[3]
            try:
                third = sys.argv[4]
                parameter_balancing_wrapper(model_name, first, second, third)
            except: parameter_balancing_wrapper(model_name, first, second)
        except: parameter_balancing_wrapper(model_name, first)
    except: parameter_balancing_wrapper(model_name)


