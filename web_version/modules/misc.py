#!/usr/bin/python
import re
import string

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
    a          = removeDoubleQuotes(sbtab_file)
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

def rerevoke_sbtab(sbtab_file):
    '''
    rerevokes the old SBtab format back to the new and current format
    '''
    rerevoked_sbtab = '!!SBtab SBtabVersion='"'1.0'"' TableType='"'Quantity'"' TableName='"'Parameter'"'\n'
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

def extract_priorpseudos(prior_file,delimiter):
    '''
    extracts the priors and pseudos of a given SBtab prior table
    '''
    pmin    = {'standard chemical potential':None,'catalytic rate constant geometric mean':None,'concentration':None,'concentration of enzyme':None,'Michaelis constant':None,'inhibitory constant':None,'activation constant':None,'chemical potential':None,'product catalytic rate constant':None,'substrate catalytic rate constant':None,'equilibrium constant':None,'forward maximal velocity':None,'reverse maximal velocity':None,'reaction affinity':None}
    pmax    = {'standard chemical potential':None,'catalytic rate constant geometric mean':None,'concentration':None,'concentration of enzyme':None,'Michaelis constant':None,'inhibitory constant':None,'activation constant':None,'chemical potential':None,'product catalytic rate constant':None,'substrate catalytic rate constant':None,'equilibrium constant':None,'forward maximal velocity':None,'reverse maximal velocity':None,'reaction affinity':None}
    priors  = {'standard chemical potential':[None,None],'catalytic rate constant geometric mean':[None,None],'concentration':[None,None],'concentration of enzyme':[None,None],'Michaelis constant':[None,None],'inhibitory constant':[None,None],'activation constant':[None,None]}
    pseudos = {'chemical potential':[None,None],'product catalytic rate constant':[None,None],'substrate catalytic rate constant':[None,None],'equilibrium constant':[None,None],'forward maximal velocity':[None,None],'reverse maximal velocity':[None,None],'reaction affinity':[None,None]}
    
    split_rows = prior_file.split('\n')
    for row in split_rows:
        split_row = row.split(delimiter)
        if split_row[0].startswith('!!') or split_row[0].startswith('%'): continue
        elif split_row[0].startswith('!'):
            for i,element in enumerate(split_row):
                if element == '!PriorMode' or element == '!PriorMedian': m_column = i
                elif element == '!PriorStd': s_column = i
                elif element == '!LowerBound': lb_column = i
                elif element == '!UpperBound': ub_column =i
        else:
            if split_row[0] in priors.keys() and split_row[0] != 'Michaelis constant':
                priors[split_row[0].lower()] = [float(split_row[m_column]),float(split_row[s_column])]
                try:
                    pmin[split_row[0].lower()] = float(split_row[lb_column])
                    pmax[split_row[0].lower()] = float(split_row[ub_column])
                except: pass
            elif split_row[0] in priors.keys():
                priors[split_row[0]] = [float(split_row[m_column]),float(split_row[s_column])]
                try:
                    pmin[split_row[0]] = float(split_row[lb_column])
                    pmax[split_row[0]] = float(split_row[ub_column])
                except: pass
            elif split_row[0] in pseudos.keys():
                pseudos[split_row[0].lower()] = [float(split_row[m_column]),float(split_row[s_column])]
                try:
                    pmin[split_row[0].lower()] = float(split_row[lb_column])
                    pmax[split_row[0].lower()] = float(split_row[ub_column])
                except: pass
            
    return priors,pseudos,pmin,pmax
    
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
    
    allowed_options = ['use_pseudos','pH','temperature','overwrite_kinetics','cell_volume','parametrisation','enzyme_prefactor','default_inhibition','default_activation','model_name','boundary_values','samples']

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
            if splitrow[o_id] == 'pH': pass
            elif not splitrow[o_id].lower() in allowed_options and splitrow[o_id] != '':
                log.append('Warning: There is an irregular option in your options file: %s'%(splitrow[o_id]))
                continue
            if splitrow[v_id] == '': continue
            voption = splitrow[o_id].lower()
            if voption == 'ph': voption = 'pH'
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
    lb_column = None
    ub_column = None
    
    split_rows = prior.split('\n')
    for row in split_rows:
        split_row = row.split(delimiter)
        if split_row[0].startswith('!!'):
            try: tabletype = re.search("TableType='([^']*)'",row).group(1)
            except:
                try: tabletype = re.search('TableType="([^"]*)"',row).group(1)
                except: tabletype = False
            if tabletype != 'QuantityInfo': validity.append('Error: The TableType of the SBtab is incorrect: "'"%s"'". It should be "'"QuantityInfo"'".'%(tabletype))
        elif split_row[0].startswith('%'): continue
        elif split_row[0].startswith('!'):
            for i,element in enumerate(split_row):
                if element == '!QuantityType': qt_column = i
                elif element == '!PriorMode' or element == '!PriorMedian': m_column = i
                elif element == '!PriorStd': s_column = i
                elif element == '!Min': lb_column = i
                elif element == '!Max': ub_column =i
        else: pass

    if qt_column == None:
        validity.append('Error: The SBtab prior table does not have the required "'"!QuantityType"'" column. Please add it to continue.')
        return prior,validity
    if m_column == None:
        validity.append('Error: The SBtab prior table does not have the required "'"!PriorMode"'" column. Please add it to continue.')
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

def removeDoubleQuotes(sbtab_file_string):
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
        

def parseReactionTable(sbtab_file,file_name,export=False):
    '''
    parses a Reaction SBtab to a stoichiometric table of the reaction formula.
    if the export parameter is set to True, the file will be written to the hard disk automatically
    '''
    from . import tablibIO
    from . import SBtab

    if sbtab_file.table_type != 'Reaction':
        print('The given TableType \"%s\" cannot be parsed. The TableType \"Reaction\" is required.'%sbtab_file.table_type)
        return False

    if not '!ReactionFormula' in sbtab_file.columns:
        print('The given provided SBtab file misses the column \"ReactionFormula\".')
        return False
    else:
        for i,c in enumerate(sbtab_file.columns):
            if c == '!ReactionFormula': rf = i
            elif c == '!ID': r = i

    react_stoich_sub  = []
    react_stoich_prod = []
    
    for row in sbtab_file.value_rows:
        reaction = row[r]
        formula  = row[rf]
        left     = formula.split('<=>')[0].lstrip().rstrip()
        right    = formula.split('<=>')[1].lstrip().rstrip()
        if '+' in left:
            subs = left.split('+')
            for sub in subs:
                sub = sub.lstrip().rstrip()
                try:
                    float(sub[0])
                    (stoich,sub) = sub.split(' ')
                    st = reaction+'\t'+stoich+'\t'+sub+'\t\t\t\n'
                    react_stoich_sub.append(st)
                except:
                    st = reaction+'\t1\t'+sub+'\t\t\t\n'
                    react_stoich_sub.append(st)
        else:
            try:
                float(left[0])
                (stoich,left) = left.split(' ')
                st = reaction+'\t'+stoich+'\t'+left+'\t\t\t\n'
                react_stoich_sub.append(st)
            except:
                st = reaction+'\t1\t'+left+'\t\t\t\n'
                react_stoich_sub.append(st)
        if '+' in right:
            prods = right.split('+')
            for prod in prods:
                prod = prod.lstrip().rstrip()
                try:
                    float(prod[0])
                    (stoich,prod) = prod.split(' ')
                    st = reaction+'\t'+stoich+'\t\t'+prod+'\t\t\n'
                    react_stoich_prod.append(st)
                except:
                    st = reaction+'\t1\t\t'+prod+'\t\n'
                    react_stoich_prod.append(st)
        else:
            try:
                float(right[0])
                (stoich,right) = right.split(' ')
                st = reaction+'\t'+stoich+'\t\t'+right+'\t\t\n'
                react_stoich_prod.append(st)
            except:
                st = reaction+'\t1\t\t'+right+'\t\t\n'
                react_stoich_prod.append(st)
        
    new_SBtab = '!!SBtab SBtabVersion="1.0" TableType="StoichiometricMatrix" TableName="%s" UniqueKey="False"\n!ReactionID\t!Stoichiometry\t!Substrate\t!Product!Location\n'%file_name[:-4]
    for sub in react_stoich_sub:
        new_SBtab += sub
    for prod in react_stoich_prod:
        new_SBtab += prod

    if export:
        parseTable = open('parseTable.tsv','w')
        parseTable.write(new_SBtab)
        parseTable.close()

    return_tab   = tablibIO.importSetNew(new_SBtab,'bla.tsv',separator='\t')
    return_SBtab = SBtab.SBtabTable(return_tab,'parseTableReaction.tsv')

    return return_SBtab
