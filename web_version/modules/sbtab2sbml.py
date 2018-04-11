"""
SBtab2SBML
==========

Python script that converts SBtab file/s to SBML.
"""
#!/usr/bin/env python
import re, libsbml
import SBtab_new
import tablibIO
import string
import random
import sys

#all allowed SBtab types
sbtab_types = ['Quantity','Event','Rule']
urns = ["obo.chebi","kegg.compound","kegg.reaction","obo.go","obo.sgd","biomodels.sbo","ec-code","kegg.orthology","uniprot"]

class ConversionError(Exception):
    '''
    Base class for errors in the SBtab conversion class.
    '''    
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return self.message

class SBtabDocument:
    '''
    SBtab document to be converted to SBML model
    '''
    def __init__(self,sbtab,filename=None,tabs=1):
        '''
        Initalizes SBtab document, checks it for SBtab count.
        If there are more than 1 SBtab file to be converted, provide a "tabs" parameter higher than 1.

        Parameters
        ----------
        sbtab : str
           SBtab file as string representation.
        filename : str
           SBtab file name.
        tabs : int
           Amount of SBtab tables in the provided file.
        '''
        self.filename  = filename

        if self.filename.endswith('tsv') or self.filename.endswith('tab') or self.filename.endswith('csv') or self.filename.endswith('.xls'): pass
        else: raise ConversionError('The given file format is not supported: '+self.filename)

        self.document = [sbtab]
        self.tabs      = tabs            
        self.unit_mM   = False
        self.unit_mpdw = False
        self.checkTabs()              #check how many SBtabs are given in the document

    def checkTabs(self):
        '''
        Checks, how many SBtab files are given by the user and saves them
        in a list, moreover stores the SBtab types in a dict linking to the SBtabs.
        '''
        self.type2sbtab = {}

        #if there are more than one SBtabs given in single files that might be comprised of several SBtabs:
        if self.tabs > 1:
            for single_document in self.document[0]:
                #check for several SBtabs in one document
                document_rows = single_document.split('\n')
                tabs_in_document = self.getAmountOfTables(document_rows)
                if tabs_in_document > 1:
                    sbtabs = self.splitDocumentInTables(document_rows)
                else: sbtabs = [document_rows]
                #generate SBtab class instance for every SBtab
                for sbtab in sbtabs:
                    sbtabtsv = self.unifySBtab(sbtab)
                    if sbtabtsv == False: continue
                    new_tablib_obj = tablibIO.importSetNew(sbtabtsv,self.filename,separator='\t')
                    single_tab = SBtab_new.SBtabTable(new_tablib_obj,self.filename)
                    if single_tab.table_type in self.type2sbtab.keys():
                        fn = random_number = str(random.randint(0,1000))
                        self.type2sbtab[single_tab.table_type+'_'+fn] = single_tab
                    else: self.type2sbtab[single_tab.table_type] = single_tab
        #elif there is only one document given, possibly consisting of several SBtabs
        else:
            #check for several SBtabs in one document
            document_rows    = self.document[0].split('\n')
            tabs_in_document = self.getAmountOfTables(document_rows)
            if tabs_in_document > 1: sbtabs = self.splitDocumentInTables(document_rows)
            else: sbtabs = [document_rows]
            #generate SBtab class instance for every SBtab
            for sbtab in sbtabs:
                as_sbtab = '\n'.join(sbtab)
                new_tablib_obj = tablibIO.importSetNew(as_sbtab,self.filename,separator='\t')
                single_tab = SBtab_new.SBtabTable(new_tablib_obj,self.filename)
                self.type2sbtab[single_tab.table_type] = single_tab

    def unifySBtab(self,sbtab):
        '''
        If we have a list of heterogeneous SBtab files, we have to unify them to one common delimiter; we choose \t arbitrarily.

        Parameters
        ----------
        sbtab : str
           SBtab file as string representation.
        '''
        new_tab = []

        for row in sbtab:
            if row.startswith('!!'): continue
            if row.startswith('!'):
                columns = row
                if '\t' in columns:
                    delimiter = '\t'
                    new_tab.append(sbtab[0])
                    new_tab.append(sbtab[1])
                    continue
                elif ';' in columns:
                    delimiter = ';'
                    new_tab.append(sbtab[0].replace(delimiter,'\t'))
                    new_tab.append(sbtab[1].replace(delimiter,'\t'))
                    continue
                elif ',' in columns:
                    delimiter = ','
                    new_tab.append(sbtab[0].replace(delimiter,'\t'))
                    new_tab.append(sbtab[1].replace(delimiter,'\t'))
                    continue
                else:
                    print 'The delimiter of one of the SBtabs could not be identified. Please check.'
            else:
                try: new_tab.append(row.replace(delimiter,'\t'))
                except: return False
            
        new_tab = '\n'.join(new_tab)

        return new_tab
            
    def getAmountOfTables(self,document_rows):
        '''
        Counts the SBtab tables that are present in the document.

        Parameters
        ----------
        document_rows : str
           Whole SBtab document as a string representation.
        '''
        counter = 0
        for row in document_rows:
            if row.startswith('!!'):
                counter += 1
        return counter

    def splitDocumentInTables(self,document_rows):
        '''
        If the document contains more than one SBtab, this function splits the document into the single SBtabs.

        Parameters
        ----------
        document_rows : str
           Whole SBtab document as a string representation.
        '''
        single_sbtab = [document_rows[0]]
        sbtab_list   = []
        for row in document_rows[1:]:
            if row.split('\t')[0] == '':
                continue
            if not row.startswith('!!'):
                    single_sbtab.append(row)
            else:
                sbtab_list.append(single_sbtab)
                single_sbtab = [row]
        sbtab_list.append(single_sbtab)
        return sbtab_list

    def makeSBML(self):
        '''
        Generates the SBML file using the provided SBtab file/s.
        '''
        # initialize new model
        self.warnings     = []
        self.new_document = libsbml.SBMLDocument()
        self.new_model    = self.new_document.createModel()
        self.new_model.setId('default_id')
        self.new_model.setName('default_name')
        self.new_document.setLevelAndVersion(2,4)
        self.reaction_list    = []
        self.species_list     = []
        self.compartment_list = []
        self.modifier_list    = []
        self.id2sbmlid        = {}
        strikes               = 1
        valid                 = True
        newSBML               = False

        while valid:
            #0st order: create compartments
            try: self.checkForCompartments()
            except:
                self.warnings.append('Error: The compartment initialisation crashed. Please check for valid compartment information.')
                break

            #1st order of bizness: due to the right modeling order of SBML, we first check for a compound SBtab
            if 'Compound' in self.type2sbtab.keys():
                try: self.compoundSBtab()
                except:
                    self.warnings.append('Warning: The provided compounds could not be initialised. Please check for valid compound information.')
            else: strikes += 1


            #2nd order of bizness: Work the Reaction SBtab (mandatory)
            if 'Reaction' in self.type2sbtab.keys():
                try:
                    self.reactionSBtab()
                except:
                    self.warnings.append('Error: The provided reaction information could not be converted. Please check for valid reaction information.')
                    break
            else: strikes += 1

            #3rd order: check, which other SBtabs are given
            for sbtab in sbtab_types:
                try:
                    self.type2sbtab[sbtab]
                    name = 'self.'+sbtab.lower()+'SBtab()'
                    eval(name)
                except:
                    pass

            #Last, but not least: generate the SBML model
            #libsbml.writeSBML(self.new_document,'New_Model.xml')
            newSBML = libsbml.writeSBMLToString(self.new_document)
            break

        if strikes < 3: return newSBML,self.warnings
        else: return False,['There was no or not sufficient model information available to build an SBML model.']

    def getWarningOnly(self):
        '''
        Returns warnings from the SBML conversion.
        '''
        return self.warnings
    
    def checkForCompartments(self):
        '''
        If there is no Compartment SBtab AND no compartments given in the other provided SBtab files, a default
        compartment needs to be set.
        '''
        self.def_comp_set = False      #has a default compartment been set?
        
        #1. check for compartment SBtab
        try:
            self.compartmentSBtab()
            return True
        except:
            pass

        #2. if there was no compartment SBtab given, check whether it is given in the other SBtabs
        try:
            sbtab = self.type2sbtab['Reaction']
            sbtab.columns_dict['!Location']
            for row in sbtab.value_rows:
                if row[sbtab.columns_dict['!Location']] != '':
                    return True
        except:
            pass

        #3. No compartment yet? Try the Compound SBtab (if present)
        try:
            sbtab = self.type2sbtab['Compound']
            sbtab.columns_dict['!Location']
            for row in sbtab.value_rows:
                if row[sbtab.columns_dict['!Location']] != '':
                    return True
        except:
            pass

        #4. Nothing yet? Then create a default compartment
        self.def_comp_set   = True
        default_compartment = self.new_model.createCompartment()
        default_compartment.setId('Default_Compartment')
        default_compartment.setName('Default_Compartment')
        default_compartment.setSize(1)
        self.compartment_list.append('Default_Compartment')
        return True

    def setAnnotation(self,element,annotation,urn,elementtype):
        '''
        Sets an annotation for a given SBML element.

        Parameters
        ----------
        element : libsbml object
           Element that needs to be annotated.
        annotation : str
           The identifier part of the annotation string.
        urn : str
           URN that links to the external web resource.
        elementtype : str
           What kind of element needs to be annotated? Model or Biological?
        '''
        element.setMetaId(element.getId()+"_meta")
        cv_term = libsbml.CVTerm()
        if elementtype == 'Model':
            cv_term.setQualifierType(0)
            cv_term.setModelQualifierType(libsbml.BQB_IS)
        else:
            cv_term.setQualifierType(1)
            cv_term.setBiologicalQualifierType(libsbml.BQB_IS)

        resource_term = "http://identifiers.org/"+urn+'/'+annotation
        cv_term.addResource(resource_term)

        return cv_term

    def compartmentSBtab(self):
        '''
        Extracts the information from the Compartment SBtab and writes it to the model.
        '''
        sbtab     = self.type2sbtab['Compartment']
        comp2size = {}

        #complement the missing compartments
        for row in sbtab.value_rows:
            if row[sbtab.columns_dict['!Compartment']] not in self.compartment_list:
                compartment = self.new_model.createCompartment()
                if '!Location:SBML:compartment:id' in sbtab.columns and row[sbtab.columns_dict['!Location:SBML:compartment:id']] != '':
                    compartment.setId(str(row[sbtab.columns_dict['!Location:SBML:compartment:id']]))
                else:
                    compartment.setId(str(row[sbtab.columns_dict['!Compartment']]))
                if '!Name' in sbtab.columns and row[sbtab.columns_dict['!Name']] != '':
                    compartment.setName(str(row[sbtab.columns_dict['!Name']]))
                else:
                    compartment.setName(str(row[sbtab.columns_dict['!Compartment']]))
                #if '!Name' in sbtab.columns and not row[sbtab.columns_dict['!Name']] == '' and not str(row[sbtab.columns_dict['!Name']]).startswith('No Name'):
                #    #if '|' in row[sbtab.columns_dict['!Name']]: compartment.setName(str(row[sbtab.columns_dict['!Name']].split('|')[0]))
                #    compartment.setName(str(row[sbtab.columns_dict['!Name']]))
                #else:
            self.compartment_list.append(row[sbtab.columns_dict['!Compartment']])

            #set the compartment sizes if given
            if '!Size' in sbtab.columns:
                for comp in self.new_model.getListOfCompartments():
                    for compsbtab in sbtab.value_rows:
                        if comp.getId() == compsbtab[sbtab.columns_dict['!Compartment']] and compsbtab[sbtab.columns_dict['!Size']] != '':
                            comp.setSize(float(compsbtab[sbtab.columns_dict['!Size']]))

            if '!SBOTerm' in sbtab.columns and row[sbtab.columns_dict['!SBOTerm']] != '':
                try: compartment.setSBOTerm(int(row[sbtab.columns_dict['!SBOTerm']][4:]))
                except: pass
            
            for column in sbtab.columns_dict.keys():
                if "Identifiers" in column:
                    annot = row[sbtab.columns_dict[column]]
                    if annot == '':
                        continue
                    for pattern in urns:
                        if pattern in column:
                            urn = pattern
                    try:
                        cv_term = self.setAnnotation(compartment,annot,urn,'Model')
                        compartment.addCVTerm(cv_term)
                    except:
                        print 'There was an annotation that I could not assign properly: ',compartment.getId(),annot #,urn
           

    def compoundSBtab(self):
        '''
        Extracts the information from the Compound SBtab and writes it to the model.
        '''
        sbtab = self.type2sbtab['Compound']

        for row in sbtab.value_rows:
            if not row[sbtab.columns_dict['!ID']] in self.species_list:
                species = self.new_model.createSpecies()
                if '!Compound:SBML:species:id' in sbtab.columns and row[sbtab.columns_dict['!Compound:SBML:species:id']] != '':
                    species.setId(str(row[sbtab.columns_dict['!Compound:SBML:species:id']]))
                    self.id2sbmlid[row[sbtab.columns_dict['!ID']]] = row[sbtab.columns_dict['!Compound:SBML:species:id']]
                else:
                    species.setId(str(row[sbtab.columns_dict['!ID']]))
                    self.id2sbmlid[row[sbtab.columns_dict['!ID']]] = None
                if '!Name' in sbtab.columns and not row[sbtab.columns_dict['!Name']] == '':
                    if '|' in row[sbtab.columns_dict['!Name']]: species.setName(str(row[sbtab.columns_dict['!Name']].split('|')[0]))
                    else: species.setName(str(row[sbtab.columns_dict['!Name']]))
                self.species_list.append(species.getId())
                #check out the speciestype if possible
                if '!SBML:speciestype:id' in sbtab.columns and row[sbtab.columns_dict['!SBML:speciestype:id']] != '':
                    species_type = self.new_model.createSpeciesType()
                    species_type.setId(str(row[sbtab.columns_dict['!SBML:speciestype:id']]))
                    species.setSpeciesType(row[sbtab.columns_dict['!SBML:speciestype:id']])
                #if compartments are given, add them
                if '!Location' in sbtab.columns and row[sbtab.columns_dict['!Location']] != '':
                    if not row[sbtab.columns_dict['!Location']] in self.compartment_list:
                        new_comp = self.new_model.createCompartment()
                        new_comp.setId(str(row[sbtab.columns_dict['!Location']]))
                        self.compartment_list.append(row[sbtab.columns_dict['!Location']])
                    species.setCompartment(row[sbtab.columns_dict['!Location']])
                elif self.def_comp_set:
                    species.setCompartment('Default_Compartment')

                if '!InitialConcentration' in sbtab.columns and row[sbtab.columns_dict['!InitialConcentration']] != '':
                    species.setInitialConcentration(float(row[sbtab.columns_dict['!InitialConcentration']]))
                elif '!InitialValue' in sbtab.columns and row[sbtab.columns_dict['!InitialValue']] != '':
                    species.setInitialConcentration(float(row[sbtab.columns_dict['!InitialValue']]))
                #DEPRECATED: Libsbml does not want this anymore!
                #is the charge at hand? if so: set
                #if sbtab.charge_column and row[sbtab.charge_column] != '':
                #    species.setCharge(int(row[sbtab.charge_column]))
                #is the species a constant and we have this information?

                if '!IsConstant' in sbtab.columns and row[sbtab.columns_dict['!IsConstant']] != '':
                    if row[sbtab.columns_dict['!IsConstant']].lower() == 'false':
                        try:
                            species.setConstant(0)
                            species.setBoundaryCondition(0)
                        except: pass
                    else:
                        try:
                            species.setConstant(1)
                            species.setBoundaryCondition(1)
                        except: pass

                if '!SBOTerm' in sbtab.columns and row[sbtab.columns_dict['!SBOTerm']] != '':
                    try: species.setSBOTerm(int(row[sbtab.columns_dict['!SBOTerm']][4:]))
                    except: pass

                if '!Unit' in sbtab.columns and self.unit_mM == False:
                    if row[sbtab.columns_dict['!Unit']] == 'mM':
                        self.makeUnitDefmM()
                        self.unit_mM = True
                    elif row[sbtab.columns_dict['!Unit']].lower().startswith('molecules'):
                        self.makeUnitDefmpdw()
                        self.unit_mpdw = True

                for column in sbtab.columns_dict.keys():
                    if "Identifiers" in column:
                        annot = row[sbtab.columns_dict[column]]
                        if annot == '': continue
                        for pattern in urns:
                            if pattern in column:
                                urn = pattern
                        try:
                            cv_term = self.setAnnotation(species,annot,urn,'Biological')
                            species.addCVTerm(cv_term)
                        except:
                            print 'There was an annotation that I could not assign properly: ',species.getId(),annot #,urn

        #species without compartments yield errors --> set them to the first available compartment
        for species in self.new_model.getListOfSpecies():
            if not species.isSetCompartment():
                species.setCompartment(self.compartment_list[0])

    def reactionSBtab(self):
        '''
        Extracts the information from the Reaction SBtab and writes it to the model.
        '''
        sbtab = self.type2sbtab['Reaction']

        #if we have the sumformulas, extract the reactants and create the species
        if '!ReactionFormula' in sbtab.columns:
            self.getReactants(sbtab)
            for reaction in self.reaction2reactants:
                try: compartment = self.reaction2compartment[reaction]
                except: compartment = False
                educts      = self.reaction2reactants[reaction][0]
                for educt in educts:
                    if educt == '': continue
                    if not educt in self.id2sbmlid.keys() and not educt in self.species_list:
                        sp = self.new_model.createSpecies()
                        sp.setId(str(educt))
                        sp.setName(str(educt))
                        sp.setInitialConcentration(1)
                        if compartment: sp.setCompartment(compartment)
                        elif self.def_comp_set: sp.setCompartment('Default_Compartment')
                        self.species_list.append(educt)
                products = self.reaction2reactants[reaction][1]
                for product in products:
                    if product == '': continue
                    if not product in self.id2sbmlid.keys() and not product in self.species_list:
                        sp = self.new_model.createSpecies()
                        sp.setId(str(product))
                        sp.setName(str(product))
                        sp.setInitialConcentration(1)
                        if compartment: sp.setCompartment(compartment)
                        elif self.def_comp_set: sp.setCompartment('Default_Compartment')
                        self.species_list.append(product)

        #if compartments are given for the reactions and these compartments are not built yet:
        if '!Location' in sbtab.columns:
            for row in sbtab.value_rows:
                if row[sbtab.columns_dict['!Location']] == '':
                    continue
                if row[sbtab.columns_dict['!Location']] not in self.compartment_list:
                    compartment = self.new_model.createCompartment()
                    compartment.setId(row[sbtab.columns_dict['!Location']])
                    compartment.setName(row[sbtab.columns_dict['!Location']])
                    compartment.setSize(1)
                    self.compartment_list.append(row[sbtab.columns_dict['!Location']])

        try:
            sbtab.columns_dict['!KineticLaw']
            self.warnings.append('Warning: Please be aware that the SBtab -> SBML conversion does not include a validation of the provided kinetic rate laws. Thus, invalid SBML code may be produced which cannot be simulated. Please check the correctness of your kinetic rate laws manually.')
        except: pass

        #creating the reactions
        for row in sbtab.value_rows:
            #if the reaction must not be included in the model: continue
            try:
                sbtab.columns_dict['!BuildReaction']
                if row[sbtab.columns_dict['!BuildReaction']] == 'False':
                    continue
            except: pass
            react = self.new_model.createReaction()
            try:
                sbtab.columns_dict['!Reaction:SBML:reaction:id']
                if row[sbtab.columns_dict['!Reaction:SBML:reaction:id']] != '':
                    react.setId(str(row[sbtab.columns_dict['!Reaction:SBML:reaction:id']]))
            except:
                react.setId(str(row[sbtab.columns_dict['!ID']]))
            if '!SBOTerm' in sbtab.columns and row[sbtab.columns_dict['!SBOTerm']] != '':
                try: react.setSBOTerm(int(row[sbtab.columns_dict['!SBOTerm']][4:]))
                except: pass

            if '!IsReversible' in sbtab.columns and row[sbtab.columns_dict['!IsReversible']] != '':
                if string.capitalize(row[sbtab.columns_dict['!IsReversible']]) == 'False':
                    try: react.setReversible(0)
                    except: pass
                elif string.capitalize(row[sbtab.columns_dict['!IsReversible']]) == 'True':
                    try: react.setReversible(1)
                    except: pass                   

            #get the reaction name
            try:
                sbtab.columns_dict['!Name']
                if not row[sbtab.columns_dict['!Name']] == '':
                    if '|' in row[sbtab.columns_dict['!Name']]: react.setName(str(row[sbtab.columns_dict['!Name']].split('|')[0]))
                    else: react.setName(str(row[sbtab.columns_dict['!Name']]))
            except: react.setName(str(row[sbtab.columns_dict['!Reaction']]))

            #if sumformula is at hand: generate reactants and products
            try:
                sbtab.columns_dict['!ReactionFormula']
                if row[sbtab.columns_dict['!ReactionFormula']] != '':
                    for educt in self.reaction2reactants[react.getId()][0]:
                        if educt == '': continue
                        reactant = react.createReactant()
                        if educt in self.id2sbmlid.keys():
                            if self.id2sbmlid[educt] != None: reactant.setSpecies(self.id2sbmlid[educt])
                            else: reactant.setSpecies(educt)
                        else: reactant.setSpecies(educt)
                        reactant.setStoichiometry(self.rrps2stoichiometry[row[sbtab.columns_dict['!ID']],educt])

                    for product in self.reaction2reactants[react.getId()][1]:
                        if product == '': continue
                        reactant = react.createProduct()
                        if product in self.id2sbmlid.keys():
                            if self.id2sbmlid[product] != None: reactant.setSpecies(self.id2sbmlid[product])
                            else: reactant.setSpecies(product)
                        else: reactant.setSpecies(product)
                        reactant.setStoichiometry(self.rrps2stoichiometry[row[sbtab.columns_dict['!ID']],product])
            except: pass
            '''
            #Uncomment, if we start using SBML Level 3, Version 1, or higher
            #if location is at hand: link reaction to the compartment
            if sbtab.columns_dict['!Location'] and row[sbtab.columns_dict['!Location']] != '':
                if not row[sbtab.columns_dict['!Location']] in self.compartment_list:
                    new_comp = self.new_model.createCompartment()
                    new_comp.setId(row[sbtab.columns_dict['!Location']])
                    react.setCompartment(new_comp)
                    self.compartment_list.append(row[sbtab.columns_dict['!Location']])
                else:
                    react.setCompartment(row[loc_column])
            '''

            #if an enzyme is given, mark it as modifier to the reaction
            try:
                sbtab.columns_dict['!Regulator']
                if row[sbtab.columns_dict['!Regulator']] != '':
                    if "|" in row[sbtab.columns_dict['!Regulator']]:
                        splits = row[sbtab.columns_dict['!Regulator']].split('|')
                        for element in splits:
                            #element = element.strip()
                            if element.startswith('+') and element[1:] in self.species_list:
                                try:
                                    mod = react.createModifier()
                                    mod.setSpecies(element[1:])
                                    mod.setSBOTerm(459)
                                    self.modifier_list.append(element)
                                except: pass
                            elif element.startswith('-') and element[1:] in self.species_list:
                                try:
                                    mod = react.createModifier()
                                    mod.setSpecies(element[1:])
                                    mod.setSBOTerm(20)
                                    self.modifier_list.append(element)
                                except: pass
                            elif element[1:] in self.species_list:
                                try:
                                    mod = react.createModifier()
                                    mod.setSpecies(element[1:])
                                    self.modifier_list.append(element)
                                    letring = str('Warning: The reaction modifier '+element+' could not be identified as either stimulator or inhibitor. Please add an SBO Term.')
                                    self.warnings.append(letring)
                                except: pass
                    else:
                        if (row[sbtab.columns_dict['!Regulator']]).startswith('+') and row[sbtab.columns_dict['!Regulator']] in self.species_list:
                            try:
                                mod = react.createModifier()
                                mod.setSpecies(row[sbtab.columns_dict['!Regulator']][1:])
                                mod.setSBOTerm(459)
                                self.modifier_list.append(row[sbtab.columns_dict['!Regulator']])
                            except: pass
                        elif (row[sbtab.columns_dict['!Regulator']]).startswith('-') and row[sbtab.columns_dict['!Regulator']] in self.species_list:
                            try:
                                mod = react.createModifier()
                                mod.setSpecies(row[sbtab.columns_dict['!Regulator']][1:])
                                mod.setSBOTerm(20)
                                self.modifier_list.append(row[sbtab.columns_dict['!Regulator']])
                            except: pass
                        elif row[sbtab.columns_dict['!Regulator']] in self.species_list:
                            try:
                                mod = react.createModifier()
                                mod.setSpecies(row[sbtab.columns_dict['!Regulator']])
                                self.modifier_list.append(row[sbtab.columns_dict['!Regulator']])
                                letring = str('Warning: The reaction modifier '+row[sbtab.columns_dict['!Regulator']]+' could not be identified as either stimulator or inhibitor. Please add an SBO Term.')
                                self.warnings.append(letring)
                            except: pass
                        self.modifier_list.append(row[sbtab.columns_dict['!Regulator']])

            except: pass

            '''
            #if metabolic regulators are given: extract them and create them
            try:
                sbtab.columns_dict['!MetabolicRegulators']
                if row[sbtab.columns_dict['!MetabolicRegulators']] != '':
                    acts,inhs = self.extractRegulators(row[sbtab.columns_dict['!MetabolicRegulators']])
                    for activator in acts:
                        if activator in self.species_list and activator not in self.modifier_list:
                            acti = react.createModifier()
                            acti.setSpecies(activator)
                            acti.setSBOTerm(459)
                            #react.addModifier(acti)
                            self.modifier_list.append(acti)
                    for inhibitor in inhs:
                        if inhibitor in self.species_list and inhibitor not in self.modifier_list:
                            inhi = react.createModifier()
                            inhi.setSpecies(inhibitor)
                            inhi.setSBOTerm(20)
                            #react.addModifier(inhi)
                            self.modifier_list.append(inhi)
            except: pass
            '''
            #set annotations if given:
            ####
            ####commented out, because by some reason this crashes the reaction (!)
            ####
            '''
            for column in sbtab.columns_dict.keys():
                if "Identifiers" in column:
                    annot = row[sbtab.columns_dict[column]]
                    if annot == '': continue
                    for pattern in urns:
                        if pattern in column:
                            urn = pattern
                    print urn
                    try:
                        cv_term = self.setAnnotation(react,annot,urn,'Biological')
                        react.addCVTerm(cv_term)
                    except:
                        print 'There was an annotation that I could not assign properly: ',react.getId(),annot #,urn

            '''            
            #since local parameters need to be entered *after* reaction creation, but *before* setting
            try:
                sbtab.columns_dict['!KineticLaw']
                if row[sbtab.columns_dict['!KineticLaw']] != '':
                    kl = react.createKineticLaw()
                    formula = row[sbtab.columns_dict['!KineticLaw']]
                    kl.setFormula(formula)
                    react.setKineticLaw(kl)
                    #for erraneous laws: remove them
                    if react.getKineticLaw().getFormula() == '':
                        react.unsetKineticLaw()
            except: pass
            
    def makeUnitDefmM(self):
        '''
        Builds unit definition; right now only in case of mM.
        '''
        ud = self.new_model.createUnitDefinition()
        ud.setId('mM')
        ud.setName('mM')

        mole = ud.createUnit()
        mole.setScale(-3)
        mole.setKind(libsbml.UNIT_KIND_MOLE)

        litre = ud.createUnit()
        litre.setExponent(-1)
        litre.setKind(libsbml.UNIT_KIND_LITRE)

    def makeUnitDefmpdw(self):
        '''
        Builds unit definition; right now only in case of molecules/gram dry weight.
        '''
        ud = self.new_model.createUnitDefinition()
        ud.setId('mmol_per_gDW_per_hr')
        ud.setName('mmol_per_gDW_per_hr')

        mole = ud.createUnit()
        mole.setScale(-3)
        mole.setExponent(1)
        mole.setKind(libsbml.UNIT_KIND_MOLE)

        litre = ud.createUnit()
        litre.setScale(0)
        litre.setExponent(-1)
        litre.setKind(libsbml.UNIT_KIND_GRAM)
       
        second = ud.createUnit()
        second.setScale(0)
        second.setExponent(-1)
        second.setMultiplier(0.000277777777777778)
        second.setKind(libsbml.UNIT_KIND_SECOND)

    def extractRegulators(self,mods):
        '''
        Extracts the regulators from the column "Regulator".

        Parameters
        ----------
        mods : str
           The modifiers of a reaction.
        '''
        activators = []
        inhibitors = []

        splits = mods.split(' ')

        for i,element in enumerate(splits):
            if element == '+': activators.append(splits[i+1])
            elif element == '-': inhibitors.append(splits[i+1])

        return activators,inhibitors

        
    def getReactants(self,sbtab):
        '''
        Extracts the reactants from the a reaction formula.

        Parameters
        ----------
        sbtab : SBtab object
           SBtab file as SBtab object.
        '''
        self.reaction2reactants   = {}
        self.rrps2stoichiometry   = {}
        self.reaction2compartment = {}
        self.specrect2compartment = {}
        educts   = []
        products = []
        
        for reaction in sbtab.value_rows:
            if '!ID' in sbtab.columns and reaction[sbtab.columns_dict['!ID']] != '': r_id = reaction[sbtab.columns_dict['!ID']]
            else: r_id   = reaction[sbtab.columns_dict['!Reaction']]
            if '!Location' in sbtab.columns: self.reaction2compartment[r_id] = reaction[sbtab.columns_dict['!Location']]
            sum_formula  = reaction[sbtab.columns_dict['!ReactionFormula']]
            #is a compartment given for the reaction? (nice, but we cannot set it (only in SBML version 3))
            if sum_formula.startswith('['):
                self.reaction2compartment[r_id] = re.search('[([^"]*)]',sum_formula).group(1)


            #check the educts
            try:
                educt_list   = re.search('([^"]*)<=>',sum_formula).group(1)
                educts       = []
                for educt in educt_list.split('+'):
                    try:
                        float(educt.lstrip().rstrip().split(' ')[0])
                        self.rrps2stoichiometry[r_id,educt.lstrip().rstrip().split(' ')[1]] = float(educt.lstrip().rstrip().split(' ')[0])
                        educts.append(educt.lstrip().rstrip().split(' ')[1])
                    except:
                        self.rrps2stoichiometry[r_id,educt.lstrip().rstrip()] = 1
                        educts.append(educt.lstrip().rstrip())
            except: pass


            #check the products
            try:
                product_list = re.search('<=>([^"]*)',sum_formula).group(1)
                products     = []
                for product in product_list.split('+'):
                    try:
                        float(product.lstrip().rstrip().split(' ')[0])
                        self.rrps2stoichiometry[r_id,product.lstrip().rstrip().split(' ')[1]] = float(product.lstrip().rstrip().split(' ')[0])
                        products.append(product.lstrip().rstrip().split(' ')[1])
                    except:
                        self.rrps2stoichiometry[r_id,product.lstrip().rstrip()] = 1
                        products.append(product.lstrip().rstrip())
            except: pass

            self.reaction2reactants[r_id] = [educts,products]

    def quantitySBtab(self):
        '''
        Extracts the information from the Quantity SBtab and writes it to the model.
        '''
        sbtab = self.type2sbtab['Quantity']

        for row in sbtab.value_rows:
            try:
                row[sbtab.columns_dict['!Description']]
                if row[sbtab.columns_dict['!Description']] == 'local parameter':
                    for reaction in self.new_model.getListOfReactions():
                        kl      = reaction.getKineticLaw()
                        formula = kl.getFormula()
                        if row[sbtab.columns_dict['!SBML:reaction:parameter:id']] in formula:
                            lp = kl.createParameter()
                            lp.setId(row[sbtab.columns_dict['!SBML:reaction:parameter:id']])
                            try: lp.setValue(float(row[sbtab.columns_dict['!Value']]))
                            except: lp.setValue(1.0)
                            try: lp.setUnits(row[sbtab.columns_dict['!Unit']])
                            except: pass
                            if '!Unit' in sbtab.columns and self.unit_mM == False:
                                if row[sbtab.columns_dict['!Unit']] == 'mM':
                                    self.makeUnitDefmM()
                                    self.unit_mM = True
                                elif row[sbtab.columns_dict['!Unit']].lower().startswith('molecules'):
                                    self.makeUnitDefmpdw()
                                    self.unit_mpdw = True
                            if '!SBOTerm' in sbtab.columns and row[sbtab.columns_dict['!SBOTerm']] != '':
                                try: lp.setSBOTerm(int(row[sbtab.columns_dict['!SBOTerm']][4:]))
                                except: pass
                else:
                    parameter = self.new_model.createParameter()
                    parameter.setId(row[sbtab.columns_dict['!SBML:reaction:parameter:id']])
                    parameter.setUnits(row[sbtab.columns_dict['!Unit']])
                    parameter.setValue(float(row[sbtab.columns_dict['!Value']]))
                    if '!SBOTerm' in sbtab.columns and row[sbtab.columns_dict['!SBOTerm']] != '':
                        try: parameter.setSBOTerm(int(row[sbtab.columns_dict['!SBOTerm']][4:]))
                        except: pass
            except:
                parameter = self.new_model.createParameter()
                parameter.setId(row[sbtab.columns_dict['!SBML:reaction:parameter:id']])
                try: parameter.setValue(float(row[sbtab.columns_dict['!Value']]))
                except: parameter.setValue(1.0)
                try: parameter.setUnits(row[sbtab.columns_dict['!Unit']])
                except: pass
                if '!SBOTerm' in sbtab.columns and row[sbtab.columns_dict['!SBOTerm']] != '':
                    try: parameter.setSBOTerm(int(row[sbtab.columns_dict['!SBOTerm']][4:]))
                    except: pass

    def eventSBtab(self):
        '''
        Extracts the information from the Event SBtab and writes it to the model.
        '''
        sbtab = self.type2sbtab['Event']

        for row in sbtab.value_rows:
            event = self.new_model.createEvent()
            event.setId(row[sbtab.columns_dict['!Event']])
            try: event.setName(row[sbtab.columns_dict['!Name']])
            except: pass
            try:
                if row[sbtab.columns_dict['!Assignments']] != '':
                    asses = row[sbtab.columns_dict['!Assignments']].split('|')
                    if len(asses) > 1:
                        for ass in asses:
                            ea  = event.createEventAssignment()
                            var = ass.split('=')[0].strip()
                            val = ass.split('=')[1].strip()
                            ea.setMath(libsbml.parseL3Formula(val))
                            ea.setVariable(var)
                    else:
                        ea  = event.createEventAssignment()
                        var = asses[0].split('=')[0].strip()
                        val = asses[0].split('=')[1].strip()
                        ea.setMath(libsbml.parseL3Formula(val))
                        ea.setVariable(var)
            except: pass            
            try:
                if row[sbtab.columns_dict['!Trigger']] != '' and row[sbtab.columns_dict['!Trigger']] != 'None':
                    trig = event.createTrigger()
                    trig.setMetaId(row[sbtab.columns_dict['!Event']]+'_meta')
                    trig.setMath(libsbml.parseL3Formula(row[sbtab.columns_dict['!Trigger']]))
            except: pass
            if '!SBOTerm' in sbtab.columns and row[sbtab.columns_dict['!SBOTerm']] != '':
                try: event.setSBOTerm(int(row[sbtab.columns_dict['!SBOTerm']][4:]))
                except: pass
            try:
                if row[sbtab.columns_dict['!Delay']] != '' and row[sbtab.columns_dict['!Delay']] != 'None':
                    dl = event.createDelay()
                    dl.setMath(libsbml.parseL3Formula(row[sbtab.columns_dict['!Delay']]))
            except: pass
            try:
                if row[sbtab.columns_dict['!UseValuesFromTriggerTime']] == 'False':
                    event.setUseValuesFromTriggerTime(0)
                else:
                    event.setUseValuesFromTriggerTime(1)
            except: pass

            for column in sbtab.columns_dict.keys():
                if "Identifiers" in column:
                    annot = row[sbtab.columns_dict[column]]
                    if annot == '': continue
                    for pattern in urns:
                        if pattern in column:
                            urn = pattern
                    try:
                        cv_term = self.setAnnotation(event,annot,urn,'Biological')
                        event.addCVTerm(cv_term)
                    except:
                        print 'There was an annotation that I could not assign properly: ',event.getId(),annot #,urn

    def ruleSBtab(self):
        '''
        Extracts the information from the Rule SBtab and writes it to the model.
        '''
        sbtab = self.type2sbtab['Rule']

        for row in sbtab.value_rows:
            if row[sbtab.columns_dict['!Name']] == 'assignmentRule':
                rule = self.new_model.createAssignmentRule()
            elif row[sbtab.columns_dict['!Name']] == 'algebraicRule':
                rule = self.new_model.createAlgebraicRule()
            elif row[sbtab.columns_dict['!Name']] == 'rateRule':
                rule = self.new_model.createRateRule()
            else: continue
            rule.setMetaId(row[sbtab.columns_dict['!Rule']]+'_meta')
            try: rule.setName(row[sbtab.columns_dict['!Name']])
            except: pass
            try: rule.setUnits(row[sbtab.columns_dict['!Unit']])
            except: pass
            try:
                if row[sbtab.columns_dict['!Formula']] != '':
                    asses = row[sbtab.columns_dict['!Formula']]
                    var = asses.split('=')[0].strip()
                    val = asses.split('=')[1].strip()
                    rule.setMath(libsbml.parseL3Formula(val))
                    rule.setVariable(var)
            except:
                pass
            for column in sbtab.columns_dict.keys():
                if "Identifiers" in column:
                    annot = row[sbtab.columns_dict[column]]
                    if annot == '': continue
                    for pattern in urns:
                        if pattern in column:
                            urn = pattern
                    try:
                        cv_term = self.setAnnotation(event,annot,urn,'Biological')
                        rule.addCVTerm(cv_term)
                    except:
                        print 'There was an annotation that I could not assign properly: ',rule.getId(),annot #,urn

if __name__ == '__main__':

    try: sys.argv[1]
    except:
        print 'You have not provided input arguments. Please start the script by also providing an SBtab file and an optional SBML output filename: >python sbtab2sbml.py SBtabfile.csv Output'
        sys.exit()
        
    file_name    = sys.argv[1]
    sbtab_file_o = open(file_name,'r')
    sbtab_file   = sbtab_file_o.read()

    try: output_name = sys.argv[2]+'.xml'
    except: output_name = file_name[:-4]+'.xml'

    Converter_class = SBtabDocument(sbtab_file,file_name)
    SBML_output     = Converter_class.makeSBML()
    new_SBML_file   = open(output_name,'w')
    new_SBML_file.write(SBML_output[0])
    new_SBML_file.close()

    print 'The SBML file has been successfully written to your working directory or chosen output path.'
