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

def table_type(sbtab):
    '''
    determines table_type of SBtab file
    '''
    for row in sbtab.split('\n'):
        if row.startswith('!!'):
            try:
                row = row.replace('"',"'")
                tabletype = re.search("TableType='([^']*)'",row).group(1)
                return tabletype
            except: pass
    return False


def validate_file_extension(file_name,file_type):
    '''
    returns Boolean to evaluate if the file has the correct extension:
    sbml => xml
    sbtab => tsv
    '''
    #check extension for sbml file
    if file_type == 'sbml' and file_name[-3:] == 'xml': return True
    elif file_type == 'sbml': return False
    else: pass

    #check extension for sbtab file
    if file_type == 'sbtab' and file_name[-3:] == 'tsv': return True
    elif file_type == 'sbtab': return False
    else: pass
    
    #if something is completely off, return False
    return False

def check_delimiter(sbtab_file):
    '''
    determine the delimiter of the tabular file
    '''
    sep = False

    for row in sbtab_file.split('\n'):
        if row.startswith('!!'): continue
        if row.startswith('!'):
            s = re.search('(.)(!)',row[1:])
            #if there is only one column, we have to define a default separator.
            #let's use a tab.
            try: sep = s.group(1)
            except: sep = '\t'

    return sep

def xml2html(sbml_file):
    '''
    generates html view out of xml file
    '''
    old_sbml = sbml_file.split('\n')
    new_sbml = '<xmp>'

    for row in old_sbml:
        new_sbml += row + '\n'

    new_sbml += '</xmp>'
        
    return new_sbml
    
def tsv2html(sbtab_file,file_name,delimiter):
    '''
    generates html view out of tsv file
    '''
    a          = remove_double_quotes(sbtab_file)
    ugly_sbtab = a.split('\n')
    nice_sbtab = '<p><h2><b>'+file_name+'</b></h2></p>'

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
            splitrow = row.split(delimiter)
        elif row.startswith('%'):
            nice_sbtab += '<tr bgcolor="#C0C0C0">'
        elif row.startswith('Parameter balancing log file'):
            nice_sbtab += '<table><tr>'
        else: nice_sbtab += '<tr>'
        
        for i,thing in enumerate(row.split(delimiter)):
            if thing.startswith('!!'): continue
            new_row = '<td>'+str(thing)+'</td>'
            nice_sbtab += new_row
        nice_sbtab += '</tr>'
    nice_sbtab += '</table>'     

    return nice_sbtab    

def revoke_sbtab(sbtab_file,delimiter):
    '''
    revokes the new SBtab format to the old SBtab format;
    this little revocation in advance is much easier than adapting the whole parameter balancing code
    to the new SBtab format.
    '''
    split_rows    = sbtab_file.split('\n')
    revoked_sbtab = ''
    for row in split_rows:
        split_row = row.split(delimiter)
        if split_row[0].startswith('!!') or split_row[0].startswith('%'): continue
        elif split_row[0].startswith('!'):
            revoked_row = []
            for element in split_row:
                revoked_element = revoke_element(element)
                revoked_row.append(revoked_element)
            revoked_sbtab += '\t'.join(revoked_row)+'\n'
        elif split_row != ['']:
            revoked_sbtab += '\t'.join(split_row)+'\n'
        
    return revoked_sbtab

def rerevoke_sbtab(sbtab_file, pseudos='no_pseudos'):
    '''
    rerevokes the old SBtab format back to the new and current format
    '''
    rerevoked_sbtab = '!!SBtab SBtabVersion='"'1.0'"' TableType='"'Quantity'"' TableName='"'Parameter_%s'"'\n' % (pseudos)
    for i,row in enumerate(sbtab_file):
        if i == 0:
            split_row = row.split('\t')            
            for elem in split_row:
                rerevoked_sbtab += rerevoke_element(elem)+'\t'
            rerevoked_sbtab = rerevoked_sbtab[:-1]+'\n'
        else:
            rerevoked_sbtab += row+'\n'

    return rerevoked_sbtab

def rerevoke_element(element):
    '''
    restore correct SBtab column headers
    '''
    rerevoke_elements = {'SBMLReactionID':'!Reaction:SBML:reaction:id','SBMLSpeciesID':'!Compound:SBML:species:id','QuantityType':'!QuantityType'}
    try: return rerevoke_elements[element]
    except: return '!'+element

def revoke_element(element):
    '''
    not only the structure, but also the column headers of the SBtab file need to be adjusted
    to the old code; this is what this little function is for
    '''
    revoke_elements = {'!Reaction:SBML:reaction:id':'SBMLReactionID','!Compound:SBML:species:id':'SBMLSpeciesID','!QuantityType':'QuantityType'}

    try: return revoke_elements[element]
    except: return element[1:]

def extract_pseudos(prior_file,delimiter):
    '''
    extracts the priors and pseudos of a given SBtab prior table
    '''
    pmin    = {'standard chemical potential':None,'catalytic rate constant geometric mean':None,'concentration':None,'concentration of enzyme':None,'Michaelis constant':None,'inhibitory constant':None,'activation constant':None,'chemical potential':None,'product catalytic rate constant':None,'substrate catalytic rate constant':None,'equilibrium constant':None,'forward maximal velocity':None,'reverse maximal velocity':None,'reaction affinity':None}
    pmax    = {'standard chemical potential':None,'catalytic rate constant geometric mean':None,'concentration':None,'concentration of enzyme':None,'Michaelis constant':None,'inhibitory constant':None,'activation constant':None,'chemical potential':None,'product catalytic rate constant':None,'substrate catalytic rate constant':None,'equilibrium constant':None,'forward maximal velocity':None,'reverse maximal velocity':None,'reaction affinity':None}
    pseudos  = {'standard chemical potential':[None]*3,'catalytic rate constant geometric mean':[None] * 3,'concentration':[None] * 3,'concentration of enzyme':[None] * 3,'Michaelis constant':[None] * 3,'inhibitory constant':[None] * 3,'activation constant':[None] * 3, 'chemical potential':[None] * 3,'product catalytic rate constant':[None] * 3,'substrate catalytic rate constant':[None] * 3,'equilibrium constant':[None] * 3,'forward maximal velocity':[None] * 3,'reverse maximal velocity':[None] * 3,'reaction affinity':[None] * 3}
    
    split_rows = prior_file.split('\n')
    for row in split_rows:
        split_row = row.split(delimiter)
        if split_row[0].startswith('!!') or split_row[0].startswith('%') or row == '': continue
        elif split_row[0].startswith('!'):
            for i,element in enumerate(split_row):
                if element == '!PriorMedian': m_column = i
                elif element == '!PriorStd': s_column = i
                elif element == '!PriorGeometricStd': pgs_column = i
                elif element == '!LowerBound': lb_column = i
                elif element == '!UpperBound': ub_column = i
                elif element == '!MathematicalType': scale_column = i
        else:
            if split_row[scale_column] == 'Additive':
                std = float(split_row[s_column])
            elif split_row[scale_column] == 'Multiplicative':
                std = float(split_row[pgs_column])
            else: std = None
            
            if split_row[0] in pseudos.keys() and split_row[0] != 'Michaelis constant':
                pseudos[split_row[0].lower()] = [float(split_row[m_column]),std]
                try:
                    pmin[split_row[0].lower()] = float(split_row[lb_column])
                    pmax[split_row[0].lower()] = float(split_row[ub_column])
                except: pass
            elif split_row[0] in pseudos.keys():
                pseudos[split_row[0]] = [float(split_row[m_column]),std]
                try:
                    pmin[split_row[0]] = float(split_row[lb_column])
                    pmax[split_row[0]] = float(split_row[ub_column])
                except: pass

    return pseudos,pmin,pmax
    
def id_checker(sbtab,sbml):
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

    for row in sbtab.split('\n'):
        splitrow = row.split('\t')
        if len(splitrow) < 3: continue
        if row.startswith('!!'): continue
        elif row.startswith('!'):
            for i,element in enumerate(splitrow):
                if element == '!Compound:SBML:species:id': s_id = i
                elif element == '!Reaction:SBML:reaction:id': r_id = i
            #continue
            if s_id == None: sbtabid2sbmlid.append('Error: The SBtab file lacks the obligatory column "'"!Compound:SBML:species:id"'" to link the parameter entries to the SBML model species.')
            if r_id == None: sbtabid2sbmlid.append('Error: The SBtab file lacks the obligatory column "'"!Reaction:SBML:reaction:id"'" to link the parameter entries to the SBML model reactions.')
        else:
            try:
                if splitrow[s_id] != '' and splitrow[s_id] not in species_ids_sbml and splitrow[s_id] != 'nan' and splitrow[s_id] != 'None':
                    sbtabid2sbmlid.append('Warning: The SBtab file holds a species ID which does not comply to any species ID in the SBML file: %s'%(splitrow[s_id]))
            except: pass
            try:
                if splitrow[r_id] != '' and splitrow[r_id] not in reaction_ids_sbml and splitrow[r_id] != 'nan' and splitrow[r_id] != 'None':
                    sbtabid2sbmlid.append('Warning: The SBtab file holds a reaction ID which does not comply to any reaction ID in the SBML file: %s'%(splitrow[r_id]))
            except: pass
            
    return sbtabid2sbmlid

def readout_config(config,delimiter):
    '''
    reads out the content of an optional config file and returns a parameter dictionary with
    many options for the balancing process
    '''
    parameter_dict = {'config':True}
    log            = []
    
    allowed_options = ['use_pseudo_values','ph','temperature','overwrite_kinetics','cell_volume','parametrisation','enzyme_prefactor','default_inhibition','default_activation','model_name','boundary_values','samples']

    for row in config.split('\n'):
        splitrow = row.split(delimiter)
        if len(splitrow) < 2: continue
        if row.startswith('!!'): continue
        elif row.startswith('!'):
            for i,element in enumerate(splitrow):
                if element == '!Option': o_id = i
                elif element == '!Value': v_id = i
            if o_id == None: log.append('Error: The config file lacks the obligatory column "'"!Option"'".')
            if v_id == None: log.append('Error: The config file lacks the obligatory column "'"!Value"'".')
        else:
            if splitrow[o_id] == 'ph': pass
            elif not splitrow[o_id].lower() in allowed_options and splitrow[o_id] != '':
                log.append('Warning: There is an irregular option in your options file: %s'%(splitrow[o_id]))
                continue
            if splitrow[v_id] == '': continue
            voption = splitrow[o_id].lower()
            vvalue  = splitrow[v_id]
            if vvalue == 'True': vvalue = True
            elif vvalue == 'False': vvalue = False
            parameter_dict[voption] = vvalue
          
    return parameter_dict,log

def cut_tabs(sbtab_file):
    '''
    cuts an SBtab document in single SBtab files
    '''
    sbtabs = {}
    current_sbtab = ''

    for row in sbtab_file.split('\n'):
        if row.startswith('!!'):
            row = row.replace('"',"'")
            if current_sbtab != '':
                sbtabs[tabletype] = current_sbtab
                tabletype = re.search("TableType='([^']*)'",row).group(1)            
                current_sbtab = row+'\n'
            else:
                tabletype = re.search("TableType='([^']*)'",row).group(1)            
                current_sbtab += row+'\n'
        else: current_sbtab += row+'\n'

    sbtabs[tabletype] = current_sbtab

    return sbtabs


def valid_prior(prior,delimiter):
    '''
    check, if the given prior file is valid and holds all required contents
    '''
    validity = []
    qts    = {'standard chemical potential':[-880,1500],'catalytic rate constant geometric mean':[10,1],'concentration':[0.1,1.5],'concentration of enzyme':[0.00001,1.5],'Michaelis constant':[0.1,1],'inhibitory constant':[0.1,1],'activation constant':[0.1,1],'chemical potential':[-880,1500],'product catalytic rate constant':[10,1.5],'substrate catalytic rate constant':[10,1.5],'equilibrium constant':[1,1.5],'forward maximal velocity':[1,2],'reverse maximal velocity':[1,2],'reaction affinity':[0,10]}
    qt_column = None
    m_column = None
    s_column = None
    logs_column = None
    lb_column = None
    ub_column = None
    
    split_rows = prior.split('\n')
    for row in split_rows:
        split_row = row.split(delimiter)
        if split_row[0].startswith('!!'):
            try: tabletype = re.search("TableType='([^']*)'",row).group(1)
            except: tabletype = False
            if tabletype != 'QuantityInfo': validity.append('Error: The TableType of the SBtab is incorrect: "'"%s"'". It should be "'"QuantityInfo"'".'%(tabletype))
        elif split_row[0].startswith('%'): continue
        elif split_row[0].startswith('!'):
            for i,element in enumerate(split_row):
                if element == '!QuantityType': qt_column = i
                elif element == '!MathematicalType': scale_column = i
                elif element == '!PriorMedian': m_column = i
                elif element == '!PriorStd': s_column = i
                elif element == '!PriorGeometricStd': pgs_column = i
                elif element == '!LowerBound': lb_column = i
                elif element == '!UpperBound': ub_column =i
        else: pass

    if qt_column == None:
        validity.append('Error: The SBtab prior table does not have the required "'"!QuantityType"'" column. Please add it to continue.')
        return prior,validity
    if m_column == None:
        validity.append('Error: The SBtab prior table does not have the required "'"!Mode"'" column. Please add it to continue.')
        return prior,validity
    if s_column == None:
        validity.append('Error: The SBtab prior table does not have the required "'"!PriorStd"'" column. Please add it to continue.')
        return prior,validity

       
    for row in split_rows:
        split_row = row.split(delimiter)
        if split_row[0].startswith('!'): pass
        elif split_row[0].startswith('%'): pass
        else:
            if split_row[0] in qts.keys():
                del qts[split_row[0]]
    
    for entry in qts.keys():
        validity.append('Warning: The SBtab prior table is missing an entry for %s. The missing value is set to the default prior distribution for this quantity.'%qts[entry])
        prior += entry+delimiter+qts[entry]+'\n'
    
    return (prior,validity)
    

def first_row(sbtab_content,delimiter):
    '''
    revokes a problem in the SBtab/tablib interface: tablib requires all rows to be
    equally long, but SBtab wants (especially for the export) only one element in
    the first row.
    '''
    splitt =  sbtab_content.split('\n')
    new_content = ''
    for i,row in enumerate(splitt):
        splitrow = row.split(delimiter)
        if i == 0: new_content += splitrow[0]+'\n'
        else: new_content += delimiter.join(splitrow)+'\n'

    return new_content

def create_filename(sbtab_name, table_type, table_name):
    '''
    creates a unique identifying name for an uploaded sbtab file to be displayed in the interface.
    '''
    if table_name != '':
        if not table_type.lower() in sbtab_name[:-4].lower(): filename = sbtab_name[:-4]+'_'+table_type+'_'+table_name
        else: filename  = sbtab_name[:-4]+'_'+table_name
    else:
        if not table_type.lower() in sbtab_name[:-4].lower(): filename = sbtab_name[:-4]+'_'+table_type
        else: filename  = sbtab_name[:-4]

    return filename

def csv2xls(sbtab_file,delimiter):
        '''
        converts sbtab file to xls file
        @sbtab_file: sbtab string
        '''
        import xlwt
        import tempfile
        
        book  = xlwt.Workbook()
        sheet = book.add_sheet('Sheet 1')

        split_sbtab_file = sbtab_file.split('\n')

        first_row = sheet.row(0)
        first_row.write(0,split_sbtab_file[0])

        for i,row in enumerate(split_sbtab_file[1:]):
            new_row   = sheet.row(i+1)
            split_row = row.split(delimiter)
            for j,element in enumerate(split_row):
                new_row.write(j,element)

        #if something is stupid and it works
        #then it is not stupid:
        book.save('simple.xls')
        fileobject = open('simple.xls','r')

        return fileobject

def remove_double_quotes(sbtab_file_string):
    '''
    remove quotes and double quotes introduced by fucking MS Excel
    '''
    try: rows = sbtab_file_string.split('\n')
    except: rows = sbtab_file_string

    sbtab = []
    for row in rows:
        n1 = row.replace('""','#')
        if n1.startswith('!!'): n2 = n1
        else: n2 = n1.replace('"','')
        new_row = n2.replace('#','"')
        sbtab.append(new_row)

    new_sbtab = '\n'.join(sbtab)

    return new_sbtab
            
def xls2csv(xls_file,filename):
        '''
        converts xls to tsv
        @xls_file: file of type xlrd
        '''
        import xlrd
        workbook = xlrd.open_workbook(filename,file_contents=xls_file)
        sheet    = workbook.sheet_by_name(workbook.sheet_names()[0])
        
        getridof = []
        csv_file = []

        for i in range(sheet.nrows):
            stringrow = str(sheet.row(i))
            notext    = string.replace(stringrow,'text:u','')
            nonumbers = string.replace(notext,'number:','')
            noopenbra = string.replace(nonumbers,'[','')
            noclosebr = string.replace(noopenbra,']','')
            if '\"!!' in noclosebr:
                noone = string.replace(noclosebr,'\"',"")
                noapostr = string.replace(noone,"\'",'\"')
            else:
                noone    = string.replace(noclosebr,"\'",'')
                noapostr = string.replace(noone,'\u201d','\"')
            nocommas  = noapostr.split(', ')
            getridof.append(nocommas)

        for row in getridof:
            new_row = ''
            for elem in row:
                if not elem == "empty:''" and not elem == 'empty:' and not elem == 'empty:""': new_row += elem+','
                else: new_row += ','
            csv_file.append(new_row.rstrip(','))

        csv_file ='\n'.join(csv_file)

        return csv_file
        
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

def is_enzyme( species ):
    return species.getSBOTerm()==14 or species.getId().startswith('enzyme') or species.getSBOTerm()==460

def get_enzyme_for_reaction(reaction, create = False):
    is_enzyme = lambda s: s.getSBOTerm()==14 or s.getId().startswith('enzyme') or s.getSBOTerm()==460
    for m in reaction.getListOfModifiers():
        s=reaction.getModel().getSpecies( m.getSpecies() )
        if is_enzyme(s):
            return s
    if create:
        e = reaction.getModel().createSpecies()
        e.setId('enzyme_'+reaction.getId())
        e.setName('enzyme_'+reaction.getId())
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
    raise Exception('No enzyme was found for reaction %s' %reaction.getId() )

def fmin_gen(f,x0,population_size=100,survivors=20,generations=20000,bounds=None,variable_is_logarithmic=None,intruders=0,use_pp=True,convenience_class=None, disp=1):
    import struct

    f = open('medians.txt','r')
    medians_no = []
    content = f.read()
    for elem in content.split(','):
        medians_no.append(float(elem))
    f.close()

    g = open('cpost.txt','r')
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
    C_post  = numpy.array(C_post_no)
    
    def local_optimize(indiv,convenience_class=None):
        better_indiv = indiv
        if convenience_class:
            better_indiv = scipy.optimize.fmin_bfgs(convenience_class.f,better_indiv,disp=0,maxiter=20)
            fval = convenience_class.f_opt(better_indiv,medians,C_post)
        else:
            #better_indiv = scipy.optimize.fmin_bfgs(f,better_indiv,disp=0,maxiter=20)
            fval = f_opt(better_indiv,medians,C_post)
        return [better_indiv,fval]
    def floatToBits(value):
        if type(value) == type(1.0):
            return (str(struct.unpack('Q', struct.pack('d', value))[0])).rjust(20,"0")
        else:
            return ",".join([floatToBits(x) for x in value.tolist()])
        
    def bitsToFloat(bits):
        if "," in bits:
            return scipy.array([bitsToFloat(x) for x in bits.split(",")])
        else:
            return struct.unpack('d', struct.pack('Q', int(bits)))[0]



    def new_individual():
        x = []
        for i in range(indiv_size):
            if variable_is_logarithmic[i]:
                logmin = scipy.log(bounds[i][0])
                logmax = scipy.log(bounds[i][1])
                x.append(scipy.exp(scipy.rand()*(logmax-logmin)+logmin))
            else:
                x.append(scipy.rand()*(bounds[i][1]-bounds[i][0])+bounds[i][0])
        return scipy.array(x)

    def bool_mate(mother,father):
        ms = floatToBits(mother)
        fs = floatToBits(father)
        l = random.sample(range(len(ms)),3)
        l.sort()
        [i1,i2,i3] = l
        cs = ms[:i1]+fs[i1:i2]+ms[i2:i3]+fs[i3:]
        child = bitsToFloat(cs)
        if child.size != mother.size or None in child or scipy.inf in child or -scipy.inf in child: raise ValueError()
        return child

    def mate(mflist):
        #return bool_mate(mflist[0],mflist[1])
        return mflist[0]+mflist[1]-mflist[2]

    def mutate(indiv):
        if not False:
            bi = floatToBits(indiv)
            number = int(math.ceil(float(len(bi))/100.))
            change_indices = random.sample([x for x in range(3,len(bi))],number)
            for ci in change_indices:
                new = bi[:ci]+str(random.choice(range(10)))+bi[(ci+1):]
                try:
                    bitsToFloat(new)
                    bi = new
                except struct.error:
                    # dont accept change
                    pass
            return scipy.absolute(bitsToFloat(bi))
        else:
            return scipy.exp(scipy.log(indiv)+wolf["mutation_factor"]*scipy.array([random.normalvariate(0,1) for x in range(indiv_size)]))

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
        bounds = [[1e-4, 1e4]]*len(x0)
    if len(bounds) != len(x0):
        raise Exception('Length of x0 and length of bounds do not fit!')

    if variable_is_logarithmic is None:
        variable_is_logarithmic = [[True]]*len(x0)
    if len(variable_is_logarithmic) != len(x0):
        raise Exception('Length of variable_is_logarithmic and x0 do not fit!')

    if use_pp:
        import pp
        # get pp servers
        servers = []
        fi = file("servers","r")
        for l in fi:
            if l.startswith("#"): continue
            servers.append(l.split(",")[1])
        fi.close()
        # read secret file
        try:
            fi = open("secret","r")
            secret = fi.readlines()[0]
            fi.close()
        except:
            import sys
            print("Please create a file calles \"secret\" with one long string in the first line.")
            sys.exit(1)
        if disp==1:
            print("starting servers",servers)
        job_server = pp.Server(ppservers=tuple(servers),secret=secret)
        print(job_server.get_active_nodes())

    population = [x0]
    indiv_size = x0.size
    quality = []
    for i in range(population_size-len(population)):
        population.append(new_individual())
    best_indiv = copy.deepcopy(x0)

    try:
        for i in range(generations):
            # rate individuals
            local_optimization_results = []
            pre_computed_qualities = len(quality)

            for j in range(pre_computed_qualities,population_size):
                if use_pp:
                    if convenience_class:
                        local_optimization_results.append( job_server.submit(local_optimize,(population[j],convenience_class),modules=("scipy","scipy.optimize","scipy.linalg","bounded"),globals=globals() ) )
                    else:
                        local_optimization_results.append( job_server.submit(local_optimize,(population[j],),modules=("scipy","scipy.optimize","scipy.linalg","bounded"),globals=globals() ) )
                else:
                    local_optimization_results.append( local_optimize(population[j]) )

            if use_pp:
                for j in range(len(local_optimization_results)):
                    local_optimization_results[j] = local_optimization_results[j]()
            for j in range(len(local_optimization_results)):
                [better_indiv,fval] = local_optimization_results[j]
                population[pre_computed_qualities+j] = better_indiv
                quality.append( fval )
            # replace None results

            for j in range(len(population)):
                while quality[j] == None or scipy.isnan(quality[j]):
                    population[j] = new_individual()
                    quality[j] = f_opt(population[j],medians,C_post)
            # sort
            sorted_quality = list(zip(quality,population))
            sorted_quality.sort(key=lambda x:x[0])
            new_population = []
            new_quality = []

            for j in range(survivors):
                f_val,indiv = sorted_quality[j]
                new_population.append(indiv)
                new_quality.append(f_val)
            population = new_population
            quality = new_quality
            # intrude
            for j in range(intruders):
                population.append(new_individual())
            # mate
            current_size = len(population)
            for j in range(current_size,population_size):
                population.append(bound(mutate(mate(random.sample(population[:current_size],3)))))
            best_indiv = copy.deepcopy(population[0])
            if disp==1:
                print("generation",str(i+1).rjust(7),"     f =",quality[0])
                print(best_indiv)

    except KeyboardInterrupt:
        if disp==1:
            print("Stopping computation")
    if use_pp:
        job_server.print_stats()
    if disp==1:
        print("generation goodbye      f =",quality[0])

    return best_indiv

def fmin_differential_evolution(f,x0,population_size=100,generations=20000,bounds=None,variable_is_logarithmic=None,crossover_factor=0.2,disp=0):
    import struct

    p = open('medians.txt','r')
    medians_no = []
    content = p.read()
    for elem in content.split(','):
        medians_no.append(float(elem))
    p.close()
    
    g = open('cpost.txt','r')
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
    C_post  = numpy.array(C_post_no)
    
    def floatToBits(value):
        if type(value) == type(1.0):
            return (str(struct.unpack('Q', struct.pack('d', value))[0])).rjust(20,"0")
        else:
            return ",".join([floatToBits(x) for x in value.tolist()])
    def bitsToFloat(bits):
        if "," not in bits:
            return struct.unpack('d', struct.pack('Q', long(bits)))[0]
        else:
            return scipy.array([bitsToFloat(x) for x in bits.split(",")])
    def new_individual():
        x = []
        for i in range(indiv_size):
            if variable_is_logarithmic[i]:
                logmin = scipy.log(bounds[i][0])
                logmax = scipy.log(bounds[i][1])
                x.append(scipy.exp(scipy.rand()*(logmax-logmin)+logmin))
            else:
                x.append(scipy.rand()*(bounds[i][1]-bounds[i][0])+bounds[i][0])
        return scipy.array(x)
    def crossover(orig,crossing_vector):
        for i in range(len(orig.tolist())):
            if random.random()<crossover_factor:
                orig[i] = crossing_vector[i]
    def bound(vector):
        for i in range(indiv_size):
            if vector[i] < bounds[i][0]:
                vector[i] = bounds[i][0]
            elif vector[i] > bounds[i][1]:
                vector[i] = bounds[i][1]

    if bounds is None:
        bounds = [[1e-4, 1e4]]*len(x0)
    if len(bounds) != len(x0):
        raise Exception('Length of x0 and length of bounds do not fit!')
    if variable_is_logarithmic is None:
        variable_is_logarithmic = [[True]]*len(x0)
    if len(variable_is_logarithmic) != len(x0):
        raise Exception('Length of variable_is_logarithmic and x0 do not fit!')

    population = [x0]
    indiv_size = x0.size
    quality = [f(x0)]
    for i in range(population_size-len(population)):
        population.append(new_individual())
        quality.append(f(population[-1]))
    best_indiv = copy.deepcopy(x0)
    try:
        for i in range(generations):
            for targetindex in range(len(population)):
                [v1,v2,v3] = random.sample(population,3)
                trial_vector = v1+v2-v3
                target_vector = population[targetindex]
                crossover(trial_vector,target_vector)
                bound(trial_vector)
                trial_quality = f(trial_vector)
                if trial_quality < quality[targetindex]:
                    population[targetindex] = trial_vector
                    quality[targetindex] = trial_quality
                
            # replace None results
            for j in range(len(population)):
                while quality[j] == None:
                    population[j] = new_individual()
                    quality[j] = f(population[j])

            best_index = quality.index(min(quality))
            if disp == 1:
                print("generation",str(i+1).rjust(7),"     f =",quality[best_index])
                print(population[best_index])
            best_indiv = copy.deepcopy(population[best_index])
    except KeyboardInterrupt:
        if disp == 1:
            print("Stopping computation")
    return best_indiv

