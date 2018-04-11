#!/usr/bin/env python
import libsbml
from . import balancer
from . import kineticizer
import os
import re
import sys
from . import misc
from . import SBtab_old

def parameter_balancing_wrapper(model_name,first=None,second=None,third=None):
    '''
    wrapper for the parameter balancing via the command line.
    receives the name of the SBML model and optional file names of SBtabs that
    need to be determined
    '''
    log_file = 'Parameter balancing log file of model %s\n\n'%(model_name)
    
    #1: open and prepare the files; then check for some rudimentary validity:
    #1.1: SBML model
    reader = libsbml.SBMLReader()
    try: sbml_file  = reader.readSBML(model_name)
    except: print('The SBML file %s could not be found.'%(sbml_name))
    try: sbml_model = sbml_file.getModel()
    except:
        print('The SBML file %s is corrupt. I quit.'%(sbml_name))
        sys.exit()
    valid_extension  = misc.validate_file_extension(model_name,'sbml')
    if not valid_extension:
        print('The SBML file %s has not the correct xml extension. I quit.'%(model_name))
        sys.exit()

    pb = balancer.ParameterBalancing(sbml_model)

    #1.2: If more arguments than the SBML file name are given, determine what SBtab
    #files we have here
    def assign_sbtab(sb,table_type,s):
        '''
        assigns the sbtab to the correct table type
        '''
        global sbtab
        global sbtab_name
        global prior
        global config
        if table_type == 'Quantity':
            sbtab      = sb
            sbtab_name = s
        elif table_type == 'QuantityInfo': prior = sb
        elif table_type == 'PbConfig': config = sb
        else: print('The provided SBtab file %s could not be assigned properly'%(table_type))

    global sbtab
    global sbtab_name
    global prior
    global config
    sbtab      = False
    sbtab_name = False
    prior      = False
    config     = False
    sbtabs     = [first,second,third]

    for s in sbtabs:
        if s != None:
            try:
                valid_extension = misc.validate_file_extension(s,'sbtab')
                if not valid_extension:
                    print('The SBtab file %s has not the correct tsv extension.'%(first))
                undetermined = open(s,'r')
                sb           = undetermined.read()
                table_type   = misc.table_type(sb)
                assign_sbtab(sb,table_type,s)
            except:
                print('The file %s could not be read and can thus not be integrated.'%(first))

    #PARAMETER FILE
    if sbtab:
        try: delimiter_sbtab = misc.check_delimiter(sb)
        except: delimiter_sbtab = '\t'
        revoked_sbtab = misc.revoke_sbtab(sbtab,delimiter_sbtab)
        use_header    = revoked_sbtab.split('\n')[0].split('\t')
        no_sbtab = False
    else: no_sbtab = True

    #PRIOR FILE
    if prior:
        try: delimiter_prior = misc.check_delimiter(prior)
        except: delimiter_prior = '\t'
        valid_prior     = misc.valid_prior(prior,delimiter_prior)
        if valid_prior != []:
            for element in valid_prior:
                log_file += str(element) +'\n'
    else:
        p = os.path.dirname(os.path.abspath(__file__))+'/files/default_files/pb_prior.tsv'
        try: prior_file = open(p,'r')
        except:
            print('The prior file (/files/default_files/pb_prior.tsv) could not be found. I quit.')
            sys.exit()
        prior = prior_file.read()
        delimiter_prior = '\t'
    (priors,pseudos,pmin,pmax) = misc.extract_priorpseudos(prior,delimiter_prior)

    #CONFIG FILE
    if config:
        delimiter_config     = misc.check_delimiter(config)
        (parameter_dict,log) = misc.readout_config(config,delimiter_config)
        if log != []:
            for element in log:
                log_file += str(element) +'\n'
    else: parameter_dict = {}

    #MAKE EMPTY SBTAB IF REQUIRED
    if no_sbtab:
        (header,rows) = pb.makeEmptySBtab(pmin,pmax,parameter_dict)
        revoked_sbtab = '\t'.join(header)+'\n'
        #for row in rows:
        #    revoked_sbtab += '\t'.join(row)+'\n'
        delimiter_sbtab = '\t'
    else:
        (header,rows)  = pb.makeSBtab(revoked_sbtab,use_header,sbtab_name,'All organisms',43,pmin,pmax,parameter_dict)
        sbtabid2sbmlid = misc.id_checker(sbtab,sbml_model)
        if sbtabid2sbmlid != []:
            for element in sbtabid2sbmlid:
                log_file.append(element)

    if not 'Temperature' in parameter_dict.keys(): parameter_dict['Temperature'] = 300
    if not 'pH' in parameter_dict.keys(): parameter_dict['pH'] = 7
    if not 'standard chemical potential' in parameter_dict.keys(): parameter_dict['standard chemical potential'] = True
    if not 'catalytic rate constant geometric mean' in parameter_dict.keys(): parameter_dict['catalytic rate constant geometric mean'] = True
    if not 'Michaelis constant' in parameter_dict.keys(): parameter_dict['Michaelis constant'] = True
    if not 'activation constant' in parameter_dict.keys(): parameter_dict['activation constant'] = True
    if not 'inhibitory constant' in parameter_dict.keys(): parameter_dict['inhibitory constant'] = True
    if not 'concentration' in parameter_dict.keys(): parameter_dict['concentration'] = True
    if not 'concentration of enzyme' in parameter_dict.keys(): parameter_dict['concentration of enzyme'] = True
    if not 'equilibrium constant' in parameter_dict.keys(): parameter_dict['equilibrium constant'] = True
    if not 'substrate catalytic rate constant' in parameter_dict.keys(): parameter_dict['substrate catalytic rate constant'] = True
    if not 'product catalytic rate constant' in parameter_dict.keys(): parameter_dict['product catalytic rate constant'] = True
    if not 'forward maximal velocity' in parameter_dict.keys(): parameter_dict['forward maximal velocity'] = True
    if not 'reverse maximal velocity' in parameter_dict.keys(): parameter_dict['reverse maximal velocity'] = True
    if not 'chemical potential' in parameter_dict.keys(): parameter_dict['chemical potential'] = True
    if not 'reaction affinity' in parameter_dict.keys(): parameter_dict['reaction affinity'] = True
    if not 'use_pseudos' in parameter_dict.keys(): parameter_dict['use_pseudos'] = False

    print('\nFiles successfully read. Start balancing.\n')
    
    #2: Parameter balancing
    if parameter_dict['use_pseudos']: filled_rows = pb.fillSBtab(rows,priors,pseudos)
    else: filled_rows = pb.fillSBtab(rows,priors)
    (mean_vector,mean_vector_inc,c_post,c_post_inc,r_matrix,shannon,(header,rows),log) = pb.makeBalancing(parameter_dict,filled_rows,revoked_sbtab,pmin,pmax)

    #3: Postprocessing
    sbtabone   = "\n".join(["\t".join([str(x) for x in row]) for row in rows])
    sbtabtwo   = pb.checkForBlanks(sbtabone)
    sbtabthree = pb.mimicMean(sbtabtwo)
    sbtabfour  = sbtabthree.split('\n')
    sbtabfive  = misc.rerevoke_sbtab(sbtabtwo.split('\n'))
    
    print('\nBalancing all model parameters! Now inserting the laws and parameters into the model\n')

    #4: inserting parameters and kinetics into SBML model
    sbtab_cl       = SBtab_old.SBtabTable(sbtabfour,sbtab_name,'QuantityType')

    try: mode = parameter_dict['param']
    except: mode = 'hal'
    try: enzyme_prefac = parameter_dict['prefac']
    except: enzyme_prefac = True
    try: def_inh = parameter_dict['default_inh']
    except: def_inh = 'complete_inh'
    try: def_act = parameter_dict['default_act']
    except: def_act = 'complete_act'
    try: overwrite = parameter_dict['overwrite']
    except: overwrite = True
    kineticizer_cs = kineticizer.kineticizer_cs(sbml_model,sbtab_cl,mode,enzyme_prefac,def_inh,def_act,True)
    
    output_name = input('\nEnter optional name for output files. If you do not enter a name, the files are named after the input file.\n')
    if output_name == '':
        try: rm = re.match('.*/(.*)',str(model_name)).group(1)[:-4]
        except: rm = str(model_name)[:-4]
        output_name = rm+'_balanced'
        
    print('\nDone... now writing output files to current working directory.\n\nGoodbye!\n\n')

    #5: Write SBtab and SBML model
    model_name = str(model_name)[:-4]
    new_sbtab_file = open(output_name+'.tsv','w')
    new_sbtab_file.write(sbtabfive)
    new_sbtab_file.close()

    sbml_code = '<?xml version="1.0" encoding="UTF-8"?>\n' + sbml_model.toSBML()
    new_sbml_model = open(output_name+'.xml','w')
    new_sbml_model.write(sbml_code)
    new_sbml_model.close()

    try:
        os.remove('cpost.txt')
        os.remove('medians.txt')
    except: pass

if __name__ == '__main__':
    
    try: model_name = sys.argv[1]
    except:
        print('\nPlease provide at least an SBML model as argument for the script, i.e. \n >python parameter_balancing.py your_sbml_file.xml \n')
        sys.exit()

    try:
        first = sys.argv[2]
        try:
            second = sys.argv[3]
            try:
                third = sys.argv[4]
                parameter_balancing_wrapper(model_name,first,second,third)
            except: parameter_balancing_wrapper(model_name,first,second)
        except: parameter_balancing_wrapper(model_name,first)
    except: parameter_balancing_wrapper(model_name)


