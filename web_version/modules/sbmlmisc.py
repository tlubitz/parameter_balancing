import libsbml
import re

'''
a collection of functions that are needed within many different modules throuhgout semanticSBML
'''

error = """%s

semanticSBML requires libSBML version 4 and the libSBML Python bindings.
Get libSBML from:

    https://sourceforge.net/projects/sbml/files/libsbml/

If libSBML is already installed, please install the Python bindings
for it as well. Otherwise, semanticSBML will not work.
"""
try:
    import libsbml
except ImportError:
    print error%'Could not import libSBML.'
else:
    try:
        if not (libsbml.LIBSBML_VERSION>=40000):
            print error%'Wrong libSBML version. Version String %s'%libsbml.LIBSBML_VERSION
    except AttributeError,e:
        print error%'No libSBML version string found.'

def setintersection(set1,set2):
    '''
    TODO replace this funciton with python native function

    this is what the python doc says:

    S.intersection(T)        S & T                    new set with elements
                                                      common to S and T
    '''
    return [e for e in set1 if e in set2]

def get_element_by_id( model, id ):
    ''' get element from libsbml model, no matter which type '''
    if model.getId()==id:
        return model
    for type in [  'Species', 'Reaction', 'Compartment', 'Rule', 'UnitDefinition', 'Event', 'Parameter']:
        ret = getattr( model, 'get'+type ) (id)
        if ret:
            return ret

def unify_reaction_scheme_delimiter(string):
    return string.replace('<->', '<=>').replace('->','<=>').replace('<-','<=>').replace(' =>',' <=>').replace('<= ','<=> ')


def parse_reaction_scheme(string):
    """ Extracts the list of reactants and products and stoichimetries from a reaction string"""
    if '(n+1)' in string:
        raise Exception('Error: Complex formations are not supported.')
    re_stoich  = re.compile( '\d+(?:\.\d+)? ' )

    string = unify_reaction_scheme_delimiter(string)
    [left,right] = string.split(' <=> ')
    educts=[]
    products=[]
    stoich_educts=[]
    stoich_products=[]
    
    for species in [ x.strip() for x in left.split(' + ')]:
        m = re_stoich.match(species)
        if m!=None:
            educts.append(species[m.end():].strip())
            stoich_educts.append(float(m.group()))
        else:
            educts.append(species.strip())
            stoich_educts.append(1)
            
    for species in right.split(' + '):
        m = re_stoich.match(species)
        if m!=None:
            products.append(species[m.end():].strip())
            stoich_products.append(float(m.group()))
        else:
            products.append(species.strip())
            stoich_products.append(1)
            
    return (educts,stoich_educts,products,stoich_products)

def make_valid_xml_id( id ):
    """ eliminate forbidden characters from designated xml id """
    if id=='':
        return ''

    id=id.strip()

    if not( id[0].isalpha() or id[0]=='_' or id[0]==':' ):
        id = '_' + id

    chars_to_replace = [(',','_'), ('"',''), ('\'',''),  ('+','plus'), ('(',''), (')',''), ('.','_'), ('-','_'), (' ', '_'), ('[',''), (']',''), (':','_') ]
    for char, replace_with in chars_to_replace:
        id = id.replace( char, replace_with )
        
    return id

def make_meta_ids( model ):
    """ 
    Assign id to metaid if metaid is not set 
    @type model: libsbml.Model
    """
    types = ['Compartments','CompartmentTypes','Constraints','Events','FunctionDefinitions','InitialAssignments','Parameters',\
             'Reactions','Rules','Species','SpeciesTypes','UnitDefinitions']
    add_types = [[model],\
                 [r.getKineticLaw() for r in model.getListOfReactions()],\
                 [kl.getListOfParameters() for kl in [r.getKineticLaw() for r in model.getListOfReactions()]]]
    for l in [getattr(model, 'getListOf'+t)() for t in types] + add_types:
        for elem in l:
            if not elem.isSetMetaId():
                elem.setMetaId(elem.getId())                            
    return model

def has_celld_annotations( document ):
	'''
	does the document contain cell designer annotations
	@rtype: bool
	'''
	return not document.toSBML().find('celldesigner') in [-1]

def remove_all_celld_annotations( document ):
    """ 
    Remove all celldesigner annotations from model 
    @type document: libsbml.Model or libsbml do
    @rtype: libsbml.Model
    """
    doc = document.toSBML()
    r1 = '(?:(?:\n?<celldesigner[^>]*?>)|(?:\n?</celldesigner[^>]*?>))*' # multiline
    r2 = '<celldesigner[^>]*?>.*?</celldesigner[^>]*?>' # single line
    doc = re.sub( r2, '', doc)
    doc = re.sub( r1, '', doc)
    d = libsbml.readSBMLFromString( '<?xml version="1.0" encoding="UTF-8"?>\n%s'%doc )
    return d

def get_modifiers(reaction):
    '''
    get iterator for the reaction modifiers
    @type reaction: libsbml.Reaction
    @rtype: Iterator
    '''
    for s in reaction.getListOfModifiers():
        yield s

def get_participants( reaction ):
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

def get_enzyme_for_reaction(reaction, create=False):
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

    
if __name__=='__main__':
    import sys
    doc = libsbml.readSBML(sys.argv[1])
    model = doc.getModel()
    m=remove_all_celld_annotations(doc)
