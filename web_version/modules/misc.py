#!/usr/bin/python
import re
import string
import libsbml
import numpy
import scipy
import scipy.linalg
import scipy.optimize
import random
import copy
import math
try:
    from . import SBtab
except:
    import SBtab


def table_type(sbtab):
    '''
    determines table_type of SBtab file
    '''
    for row in sbtab.split('\n'):
        if row.startswith('!!'):
            try:
                row = row.replace('"', "'")
                tabletype = re.search("TableType='([^']*)'", row).group(1)
                return tabletype
            except: pass
    return False


def validate_file_extension(file_name, file_type):
    '''
    returns Boolean to evaluate if the file has the correct extension:
    sbml => xml
    sbtab => tsv
    '''
    # check extension for sbml file
    if file_type == 'sbml' and file_name[-3:] == 'xml': return True
    elif file_type == 'sbml': return False
    else: pass

    # check extension for sbtab file
    if file_type == 'sbtab' and file_name[-3:] == 'tsv': return True
    elif file_type == 'sbtab': return False
    else: pass

    # if something is completely off, return False
    return False


def check_delimiter(sbtab_file):
    '''
    determine the delimiter of the tabular file
    '''
    sep = False

    for row in sbtab_file.split('\n'):
        if row.startswith('!!'): continue
        if row.startswith('!'):
            s = re.search('(.)(!)', row[1:])
            # if there is only 1 column, we have to define a default separator
            # let's use a tab.
            try: sep = s.group(1)
            except: sep = '\t'

    return sep


def valid_prior(sbtab_prior):
    '''
    if the given SBtab file is a prior for parameter balancing, it needs to be
    checked thorougly for the validity of several features
    '''
    validity = []

    # check table type
    if sbtab_prior.table_type != 'Quantity':
        validity.append('Error: The TableType of the prior file is not '\
                        'correct: %s. '\
                        'It should be Quantity' % sbtab_prior.table_type)

    # check for required columns
    required_columns = ['!QuantityType', '!Unit', '!MathematicalType',
                        '!PriorMedian', '!PriorStd', '!PriorGeometricStd',
                        '!DataStd', '!Dependence',
                        '!UseAsPriorInformation', '!MatrixInfo']
    for column in required_columns:
        if column not in sbtab_prior.columns_dict:
            validity.append('Error: The crucial column %s is missing in'\
                            ' the prior file.' % column)

    # check for required row entries
    required_rows = ['standard chemical potential',
                     'catalytic rate constant geometric mean',
                     'concentration', 'concentration of enzyme',
                     'Michaelis constant', 'inhibitory constant',
                     'activation constant', 'chemical potential',
                     'product catalytic rate constant',
                     'substrate catalytic rate constant',
                     'equilibrium constant', 'forward maximal velocity',
                     'reverse maximal velocity', 'reaction affinity']
    for row in sbtab_prior.value_rows:
        try: required_rows.remove(row[sbtab_prior.columns_dict['!QuantityType']])
        except: pass

    for row in required_rows:
        validity.append('Error: The prior file is missing an entry for th'\
                        'crucial value %s.' % row)

    return validity


def extract_pseudos_priors(sbtab_prior):
    '''
    extracts the priors and pseudos of a given SBtab prior table
    '''
    pseudo_list = ['chemical potential', 'product catalytic rate constant',
                   'substrate catalytic rate constant',
                   'equilibrium constant', 'forward maximal velocity',
                   'reverse maximal velocity', 'reaction affinity']
    pmin = {}
    pmax = {}
    pseudos = {}
    priors = {}
    
    for row in sbtab_prior.value_rows:
        pmin[row[sbtab_prior.columns_dict['!QuantityType']]] = float(row[sbtab_prior.columns_dict['!LowerBound']])
        pmax[row[sbtab_prior.columns_dict['!QuantityType']]] = float(row[sbtab_prior.columns_dict['!UpperBound']])
        if row[sbtab_prior.columns_dict['!MathematicalType']] == 'Additive':
            std = row[sbtab_prior.columns_dict['!PriorStd']]
        else:
            std = row[sbtab_prior.columns_dict['!PriorGeometricStd']]
        median = row[sbtab_prior.columns_dict['!PriorMedian']]

        if row[sbtab_prior.columns_dict['!QuantityType']] in pseudo_list:
            pseudos[row[sbtab_prior.columns_dict['!QuantityType']]] = [float(median),
                                                                       float(std)]
        else:
            priors[row[sbtab_prior.columns_dict['!QuantityType']]] = [float(median),
                                                                      float(std)]
       
    return pseudos, priors, pmin, pmax


def id_checker(sbtab, sbml):
    '''
    this function checks, whether all the entries of the SBML ID columns of
    the SBtab file can also be found in the SBML file. If not, these are
    omitted during the balancing. But there should be a warning to raise user
    awareness.
    '''
    sbtabid2sbmlid = []
    reaction_ids_sbml = []
    species_ids_sbml = []
    s_id = None
    r_id = None

    for reaction in sbml.getListOfReactions():
        reaction_ids_sbml.append(reaction.getId())
    for species in sbml.getListOfSpecies():
        species_ids_sbml.append(species.getId())

    for row in sbtab.split('\n'):
        splitrow = row.split('\t')
        if len(splitrow) < 3: continue
        if row.startswith('!!'): continue
        elif row.startswith('!'):
            for i, element in enumerate(splitrow):
                if element == '!Compound:SBML:species:id': s_id = i
                elif element == '!Reaction:SBML:reaction:id': r_id = i
            if s_id is None:
                sbtabid2sbmlid.append('''Error: The SBtab file lacks the obliga
                tory column "'"!Compound:SBML:species:id"'" to link the paramet
                er entries to the SBML model species.''')
            if r_id is None:
                sbtabid2sbmlid.append('''Error: The SBtab file lacks the obliga
                tory column "'"!Reaction:SBML:reaction:id"'" to link the parame
                ter entries to the SBML model reactions.''')
        else:
            try:
                if splitrow[s_id] != '' \
                   and splitrow[s_id] not in species_ids_sbml \
                   and splitrow[s_id] != 'nan' and splitrow[s_id] != 'None':
                    sbtabid2sbmlid.append('''Warning: The SBtab file holds a sp
                    ecies ID which does not comply to any species ID in the SBM
                    L file: %s''' % (splitrow[s_id]))
            except: pass
            try:
                if splitrow[r_id] != '' \
                   and splitrow[r_id] not in reaction_ids_sbml \
                   and splitrow[r_id] != 'nan' and splitrow[r_id] != 'None':
                    sbtabid2sbmlid.append('''Warning: The SBtab file holds a re
                    action ID which does not comply to any reaction ID in the S
                    BML file: %s''' % (splitrow[r_id]))
            except: pass

    return sbtabid2sbmlid


def readout_config(sbtab_options):
    '''
    reads out the content of an optional config file and returns a parameter
    dictionary with many options for the balancing process
    '''
    parameter_dict = {'config': True}
    log = []
    allowed_options = ['use_pseudo_values', 'ph', 'temperature',
                       'overwrite_kinetics', 'cell_volume', 'parametrisation',
                       'enzyme_prefactor', 'default_inhibition',
                       'default_activation', 'model_name', 'boundary_values',
                       'samples']


    if '!Option' not in sbtab_options.columns_dict:
        log.append('Error: The crucial option column is missing from the'\
                   'options file')
    if '!Value' not in sbtab_options.columns_dict:
            log.append('Error: The crucial value column is missing from the'\
                       'options file')

    for row in sbtab_options.value_rows:
        if row[sbtab_options.columns_dict['!Option']] not in allowed_options:
            log.append('There was an irregular option in the options file:'\
                       '%s' % row[sbtab_options.columns_dict['!Option']])
        else:
            if row[sbtab_options.columns_dict['!Value']] == '':
                log.append('There is no value set for option:'\
                           '%s' % row[sbtab_options.columns_dict['!Option']])
    
        parameter_dict[row[sbtab_options.columns_dict['!Option']]] = row[sbtab_options.columns_dict['!Value']]

    return parameter_dict, log


def get_modifiers(reaction):
    '''
    get iterator for the reaction modifiers
    @type reaction: libsbml.Reaction
    @rtype: Iterator
    '''
    for s in reaction.getListOfModifiers():
        yield s


def get_participants(reaction):
    """
    get iterator for the reaction participants (reactants + products)
    @type reaction: libsbml.Reaction
    @rtype: Iterator
    """
    for s in reaction.getListOfReactants():
        yield s
    for s in reaction.getListOfProducts():
        yield s


def is_enzyme(species):
    return species.getSBOTerm() == 14 or species.getId().startswith('enzyme') \
        or species.getSBOTerm() == 460


def get_enzyme_for_reaction(reaction, create=False):
    is_enzyme = lambda s: s.getSBOTerm() == 14 \
                or s.getId().startswith('enzyme') or s.getSBOTerm() == 460
    for m in reaction.getListOfModifiers():
        s = reaction.getModel().getSpecies(m.getSpecies())
        if is_enzyme(s):
            return s
    if create:
        e = reaction.getModel().createSpecies()
        e.setId('enzyme_' + reaction.getId())
        e.setName('enzyme_' + reaction.getId())
        try:
            comp = reaction.getModel().getSpecies(reaction.getReactant(0).getSpecies()).getCompartment()
        except:
            comp = reaction.getModel().getSpecies(reaction.getProduct(0).getSpecies()).getCompartment()
        if not comp:
            comp = reaction.getModel().getCompartment(0).getId()
        e.setCompartment(comp)
        e.setSBOTerm(460)
        mod = reaction.createModifier()
        mod.setSpecies(e.getId())
        return e
    raise Exception('No enzyme was found for reaction %s' % reaction.getId())


def fmin_gen(f, x0, population_size=100, survivors=20, generations=20000,
             bounds=None, variable_is_logarithmic=None, intruders=0,
             use_pp=True, convenience_class=None, disp=1):
    import struct

    f = open('medians.txt', 'r')
    medians_no = []
    content = f.read()
    for elem in content.split(','):
        medians_no.append(float(elem))
    f.close()

    g = open('cpost.txt', 'r')
    content2 = g.read()
    C_post_no = []
    for line in content2.split('\n'):
        if line != '':
            linecontent = line.split(',')
            single_line = []
            for elem in linecontent:
                single_line.append(float(elem))
            C_post_no.append(single_line)

    medians = numpy.array(medians_no)
    C_post = numpy.array(C_post_no)

    def local_optimize(indiv, convenience_class=None):
        better_indiv = indiv
        if convenience_class:
            better_indiv = scipy.optimize.fmin_bfgs(convenience_class.f,
                                                    better_indiv, disp=0,
                                                    maxiter=20)
            fval = convenience_class.f_opt(better_indiv, medians, C_post)
        else:
            fval = f_opt(better_indiv, medians, C_post)
        return [better_indiv, fval]

    def float_to_bits(value):
        if type(value) == float:
            return (str(struct.unpack('Q',
                                      struct.pack('d',
                                                  value))[0])).rjust(20, "0")
        else:
            return ",".join([float_to_bits(x) for x in value.tolist()])

    def bits_to_float(bits):
        if "," in bits:
            return scipy.array([bits_to_float(x) for x in bits.split(",")])
        else:
            return struct.unpack('d', struct.pack('Q', int(bits)))[0]

    def new_individual():
        x = []
        for i in range(indiv_size):
            if variable_is_logarithmic[i]:
                logmin = scipy.log(bounds[i][0])
                logmax = scipy.log(bounds[i][1])
                x.append(scipy.exp(scipy.rand() * (logmax - logmin) + logmin))
            else:
                x.append(scipy.rand() * (bounds[i][1] - bounds[i][0]) +
                         bounds[i][0])
        return scipy.array(x)

    def bool_mate(mother, father):
        ms = float_to_bits(mother)
        fs = float_to_bits(father)
        l = random.sample(range(len(ms)), 3)
        l.sort()
        [i1, i2, i3] = l
        cs = ms[:i1] + fs[i1:i2] + ms[i2:i3] + fs[i3:]
        child = bits_to_float(cs)
        if child.size != mother.size or None in child or scipy.inf in child \
           or -scipy.inf in child: raise ValueError()
        return child

    def mate(mflist):
        return mflist[0] + mflist[1] - mflist[2]

    def mutate(indiv):
        if not False:
            bi = float_to_bits(indiv)
            number = int(math.ceil(float(len(bi)) / 100.))
            change_indices = random.sample([x for x in range(3, len(bi))],
                                           number)
            for ci in change_indices:
                new = bi[:ci] + str(random.choice(range(10))) + bi[(ci + 1):]
                try:
                    bits_to_float(new)
                    bi = new
                except struct.error:
                    # dont accept change
                    pass
            return scipy.absolute(bits_to_float(bi))
        else:
            return scipy.exp(scipy.log(indiv) + wolf["mutation_factor"] *
                             scipy.array([random.normalvariate(0, 1) for x in range(indiv_size)]))

    def bound(vector):
        global correct_vector

        if len(vector) == indiv_size:
            correct_vector = vector

        if len(vector) != indiv_size:
            print('im doing it right now!')
            vector = correct_vector

        for i in range(indiv_size):
            if vector[i] < bounds[i][0]:
                vector[i] = bounds[i][0]
                correct_vector = vector
            elif vector[i] > bounds[i][1]:
                vector[i] = bounds[i][1]
                correct_vector = vector
        return vector

    if bounds is None:
        bounds = [[1e-4, 1e4]] * len(x0)
    if len(bounds) != len(x0):
        raise Exception('Length of x0 and length of bounds do not fit!')
    if variable_is_logarithmic is None:
        variable_is_logarithmic = [[True]] * len(x0)
    if len(variable_is_logarithmic) != len(x0):
        raise Exception('Length of variable_is_logarithmic and x0 do not fit!')

    if use_pp:
        import pp
        # get pp servers
        servers = []
        fi = file("servers", "r")
        for l in fi:
            if l.startswith("#"): continue
            servers.append(l.split(",")[1])
        fi.close()
        # read secret file
        try:
            fi = open("secret", "r")
            secret = fi.readlines()[0]
            fi.close()
        except:
            import sys
            print('''Please create a file called \"secret\" with one long
            string in the first line.''')
            sys.exit(1)
        if disp == 1:
            print("starting servers", servers)
        job_server = pp.Server(ppservers=tuple(servers), secret=secret)
        print(job_server.get_active_nodes())

    population = [x0]
    indiv_size = x0.size
    quality = []
    for i in range(population_size - len(population)):
        population.append(new_individual())
    best_indiv = copy.deepcopy(x0)

    try:
        for i in range(generations):
            # rate individuals
            local_optimization_results = []
            pre_computed_qualities = len(quality)

            for j in range(pre_computed_qualities, population_size):
                if use_pp:
                    if convenience_class:
                        local_optimization_results.append(job_server.submit(local_optimize,
                                                                            (population[j],
                                                                             convenience_class),
                                                                            modules=("scipy",
                                                                                     "scipy.optimize",
                                                                                     "scipy.linalg",
                                                                                     "bounded"),
                                                                            globals=globals()))
                    else:
                        local_optimization_results.append(job_server.submit(local_optimize,
                                                                            (population[j],),
                                                                            modules=("scipy",
                                                                                     "scipy.optimize",
                                                                                     "scipy.linalg",
                                                                                     "bounded"),
                                                                            globals=globals()))
                else:
                    local_optimization_results.append(local_optimize(population[j]))

            if use_pp:
                for j in range(len(local_optimization_results)):
                    local_optimization_results[j] = local_optimization_results[j]()
            for j in range(len(local_optimization_results)):
                [better_indiv, fval] = local_optimization_results[j]
                population[pre_computed_qualities + j] = better_indiv
                quality.append(fval)

            # replace None results
            for j in range(len(population)):
                while quality[j] is None or scipy.isnan(quality[j]):
                    population[j] = new_individual()
                    quality[j] = f_opt(population[j], medians, C_post)

            # sort
            sorted_quality = list(zip(quality, population))
            sorted_quality.sort(key=lambda x: x[0])
            new_population = []
            new_quality = []
            for j in range(survivors):
                f_val, indiv = sorted_quality[j]
                new_population.append(indiv)
                new_quality.append(f_val)
            population = new_population
            quality = new_quality

            # intrude
            for j in range(intruders):
                population.append(new_individual())

            # mate
            current_size = len(population)
            for j in range(current_size, population_size):
                population.append(bound(mutate(mate(random.sample(population[:current_size],
                                                                  3)))))
            best_indiv = copy.deepcopy(population[0])
            if disp == 1:
                print("generation", str(i + 1).rjust(7), "     f =",
                      quality[0])
                print(best_indiv)
    except KeyboardInterrupt:
        if disp == 1:
            print("Stopping computation")

    if use_pp:
        job_server.print_stats()
    if disp == 1:
        print("generation goodbye      f =", quality[0])

    return best_indiv


def fmin_differential_evolution(f, x0, population_size=100, generations=20000,
                                bounds=None, variable_is_logarithmic=None,
                                crossover_factor=0.2, disp=0):
    import struct

    p = open('medians.txt', 'r')
    medians_no = []
    content = p.read()
    for elem in content.split(','):
        medians_no.append(float(elem))
    p.close()

    g = open('cpost.txt', 'r')
    content2 = g.read()
    C_post_no = []
    for line in content2.split('\n'):
        if line != '':
            linecontent = line.split(',')
            single_line = []
            for elem in linecontent:
                single_line.append(float(elem))
            C_post_no.append(single_line)

    medians = numpy.array(medians_no)
    C_post = numpy.array(C_post_no)

    def float_to_bits(value):
        if type(value) == float:
            return (str(struct.unpack('Q', struct.pack('d',
                                                       value))[0])).rjust(20,
                                                                          "0")
        else:
            return ",".join([float_to_bits(x) for x in value.tolist()])

    def bits_to_float(bits):
        if "," not in bits:
            return struct.unpack('d', struct.pack('Q', long(bits)))[0]
        else:
            return scipy.array([bits_to_float(x) for x in bits.split(",")])

    def new_individual():
        x = []
        for i in range(indiv_size):
            if variable_is_logarithmic[i]:
                logmin = scipy.log(bounds[i][0])
                logmax = scipy.log(bounds[i][1])
                x.append(scipy.exp(scipy.rand() * (logmax - logmin) + logmin))
            else:
                x.append(scipy.rand() * (bounds[i][1] - bounds[i][0]) +
                         bounds[i][0])
        return scipy.array(x)

    def crossover(orig, crossing_vector):
        for i in range(len(orig.tolist())):
            if random.random() < crossover_factor:
                orig[i] = crossing_vector[i]

    def bound(vector):
        for i in range(indiv_size):
            if vector[i] < bounds[i][0]:
                vector[i] = bounds[i][0]
            elif vector[i] > bounds[i][1]:
                vector[i] = bounds[i][1]

    if bounds is None:
        bounds = [[1e-4, 1e4]] * len(x0)
    if len(bounds) != len(x0):
        raise Exception('Length of x0 and length of bounds do not fit!')
    if variable_is_logarithmic is None:
        variable_is_logarithmic = [[True]] * len(x0)
    if len(variable_is_logarithmic) != len(x0):
        raise Exception('Length of variable_is_logarithmic and x0 do not fit!')

    population = [x0]
    indiv_size = x0.size
    quality = [f(x0)]

    for i in range(population_size - len(population)):
        population.append(new_individual())
        quality.append(f(population[-1]))
    best_indiv = copy.deepcopy(x0)

    try:
        for i in range(generations):
            for targetindex in range(len(population)):
                [v1, v2, v3] = random.sample(population, 3)
                trial_vector = v1 + v2 - v3
                target_vector = population[targetindex]
                crossover(trial_vector, target_vector)
                bound(trial_vector)
                trial_quality = f(trial_vector)
                if trial_quality < quality[targetindex]:
                    population[targetindex] = trial_vector
                    quality[targetindex] = trial_quality

            # replace None results
            for j in range(len(population)):
                while quality[j] is None:
                    population[j] = new_individual()
                    quality[j] = f(population[j])

            best_index = quality.index(min(quality))
            if disp == 1:
                print("generation", str(i + 1).rjust(7), "     f =",
                      quality[best_index])
                print(population[best_index])
            best_indiv = copy.deepcopy(population[best_index])

    except KeyboardInterrupt:
        if disp == 1:
            print("Stopping computation")

    return best_indiv


def xml2html(sbml_file):
    '''
    generates html view out of xml file
    '''
    old_sbml = str(sbml_file).split('\\n')
    new_sbml = '<xmp>'

    for row in old_sbml:
        new_sbml += row + '\n'

    new_sbml += '</xmp>'
        
    return new_sbml


def id_checker(sbtab, sbml):
    '''
    this function checks, whether all the entries of the SBML ID columns of the SBtab file can also be
    found in the SBML file. If not, these are omitted during the balancing. But there should be a warning
    to raise user awareness.
    '''
    sbtabid2sbmlid = []

    reaction_ids_sbml = []
    species_ids_sbml  = []

    s_id = None
    r_id = None

    for reaction in sbml.getListOfReactions():
        reaction_ids_sbml.append(reaction.getId())
    for species in sbml.getListOfSpecies():
        species_ids_sbml.append(species.getId())

    for row in sbtab.value_rows:
        if len(row) < 3: continue
        try:
            s_id = sbtab.columns_dict['!Compound:SBML:species:id']
            r_id = sbtab.columns_dict['!Reaction:SBML:reaction:id']
        except:
            sbtabid2sbmlid.append('Error: The SBtab file lacks either of the obligatory columns "'"!Compound:SBML:species:id"'" or "'"!Reaction:SBML:reaction:id"'" to link the parameter entries to the SBML model species.')
        try:
            if row[s_id] != '' and row[s_id] not in species_ids_sbml and row[s_id] != 'nan' and row[s_id] != 'None':
                sbtabid2sbmlid.append('Warning: The SBtab file holds a species ID which does not comply to any species ID in the SBML file: %s'%(row[s_id]))
        except: pass
        try:
            if row[r_id] != '' and row[r_id] not in reaction_ids_sbml and row[r_id] != 'nan' and row[r_id] != 'None':
                sbtabid2sbmlid.append('Warning: The SBtab file holds a reaction ID which does not comply to any reaction ID in the SBML file: %s'%(row[r_id]))
        except: pass
            
    return sbtabid2sbmlid


def cut_tabs(sbtab_strings):
    '''
    cuts an SBtab document in single SBtab files
    '''
    sbtabs = []
    sbtab_string = ''
    counter = 1

    for row in sbtab_strings.split('\n'):
        if row.startswith('!!'):
            if sbtab_string == '':
                sbtab_string = row + '\n'
                continue
            else:
                try:
                    sbtabs.append(sbtab_string)
                    sbtab_string = row + '\n'
                    counter += 1
                except:
                    print('Warning: Could not write SBtab %s' % counter)
                    counter += 1
        else:
            sbtab_string += row + '\n'
    sbtabs.append(sbtab_string) 
                    
    return sbtabs


    
    '''
    # THIS HERE DOESN'T HELP US WITHOUT THE SBTAB DOCUMENT CLASS
    sbtabs = []
    sbtab_string = ''
    counter = 1

    for row in sbtab_strings.split('\n'):
        print(row)
        if row.startswith('!!'):
            if sbtab_string == '':
                sbtab_string = row + '\n'
                continue
            else:
                try:
                    sbtab = SBtab.SBtabTable(sbtab_string, 'sbtab_%s.tsv' % counter)
                    sbtabs.append(sbtab)
                    sbtab_string = row + '\n'
                    print('Tis written')
                    counter += 1
                except:
                    print('Warning: Could not write SBtab %s' % counter)
                    counter += 1
        else:
            sbtab_string += row + '\n'

    sbtab = SBtab.SBtabTable(sbtab_string, 'sbtab_%s.tsv' % counter)
    sbtabs.append(sbtab) 
                    
    return sbtabs
    '''

def tsv_to_html(sbtab):
    '''
    generates html view out of tsv file
    '''
    ugly_sbtab = sbtab.return_table_string().split('\n')
    nice_sbtab = '<p><h2><b>'+sbtab.filename+'</b></h2></p>'

    first = True
    for row in ugly_sbtab:
        if row.startswith('!!') and first:
            nice_sbtab += '<a style="background-color:#00BFFF">'+row+'</a><br>'
            nice_sbtab += '<table>'
            first = False
        elif row.startswith('!!'):
            nice_sbtab += '</table>'
            nice_sbtab += '<a style="background-color:#00BFFF">'+row+'</a><br>'
            nice_sbtab += '<table>'
        elif row.startswith('!'):
            nice_sbtab += '<tr bgcolor="#87CEFA">'
            splitrow = row.split(sbtab.delimiter)
        elif row.startswith('%'):
            nice_sbtab += '<tr bgcolor="#C0C0C0">'
        elif row.startswith('Parameter balancing log file'):
            nice_sbtab += '<table><tr>'
        else: nice_sbtab += '<tr>'
        
        for i,thing in enumerate(row.split(sbtab.delimiter)):
            if thing.startswith('!!'): continue
            new_row = '<td>'+str(thing)+'</td>'
            nice_sbtab += new_row
        nice_sbtab += '</tr>'
    nice_sbtab += '</table>'     

    return nice_sbtab  
