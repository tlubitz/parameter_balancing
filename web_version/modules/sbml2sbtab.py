"""
SBML2SBtab Converter
====================

Python script that converts an SBML file to SBtab file/s.

See specification for further information.
"""
#!/usr/bin/env python
import re, libsbml, numpy
import SBtab
import sys

allowed_sbtabs = ['Compartment','Compound','Reaction','Rule','Quantity','Event']

class ConversionError(Exception):
    '''
    Base class for errors in the SBtab conversion class.
    '''
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return self.message

class SBMLDocument:
    '''
    SBML model to be converted to SBtab file/s
    '''
    def __init__(self,sbml_model,filename):
        '''
        Initalizes SBtab document, checks it for SBtabs

        Parameters
        ----------
        sbml_model : libsbml object
            SBML model as libsbml object.
        filename : str
            Filename with extension.
        '''
        #reader      = libsbml.SBMLReader()
        #sbml_modelx = reader.readSBML(sbml_model)
        #print sbml_modelx.getModel()

        self.model    = sbml_model
        self.filename = filename
        if not self.filename.endswith('.xml') and not filename.endswith('.sbml'): 
            raise ConversionError('The given file format is not supported: '+self.filename)

        #for testing purposes:
        '''
        sbtabs = self.makeSBtabs()
        for sbtab in sbtabs:
            print sbtab
            print '\n\n\n'
        '''

    def makeSBtabs(self):
        '''
        Generates the SBtab files.
        '''
        self.warnings = []
        sbtabs        = []

        #self.testForInconvertibles()

        for sbtab_type in allowed_sbtabs:
            try:
                function_name = 'self.'+sbtab_type.lower()+'SBtab()'
                new_sbtab = eval(function_name)
                if new_sbtab != False: sbtabs.append(new_sbtab)
            except:
                pass

        sbtabs = self.getRidOfNone(sbtabs)
        sbtabs = self.getRidOfEmptyColumns(sbtabs)

        sbtab_objects = []
        
        for sbtab in sbtabs:
            so = SBtab.SBtabTable(sbtab[0], self.filename[:-4]+'_%s.tsv' % sbtab[1])
            sbtab_objects.append(so)
        
        return (sbtab_objects, self.warnings)

    def testForInconvertibles(self):
        '''
        Some SBML entities cannot be converted to SBML. Thus, we add warnings to tell the user about omitting them.
        '''
        try:
            rules = self.model.getListOfRules()
            if len(rules)>0:
                self.warnings.append('The SBML model contains rules. These cannot be translated to the SBtab files yet.')
        except: pass
        
    def getRidOfNone(self,sbtabs):
        '''
        Removes empty SBtabs (if no values for an SBtab were provided by the SBML).

        Parameters
        ----------
        sbtabs : list
           List of single SBtab files.
        '''
        new_tabs = []
        for element in sbtabs:
            if element != None:
                new_tabs.append(element)
        return new_tabs

    def getRidOfEmptyColumns(self,sbtabs):
        '''
        Removes empty columns.

        Parameters
        ----------
        sbtabs : list
           List of single SBtab files.        
        '''
        sbtab2empty_columns = {}

        for i,sbtab in enumerate(sbtabs):
            tab   = sbtab[0]
            sptab = tab.split('\n')[1]
            empty_columns = list(range(0,len(sptab.split('\t'))))
            for k,element in enumerate(empty_columns): empty_columns[k] = str(element)
            for row in sbtab[0].split('\n')[2:]:
                for j,element in enumerate(row.split('\t')):
                    if element != '':
                        try: empty_columns.remove(str(j))
                        except: pass
            for k,element in enumerate(empty_columns): empty_columns[k] = int(element)
            sbtab2empty_columns[i] = empty_columns

        new_tabs = []
        
        for i,sbtab in enumerate(sbtabs):
            eliminate_columns = sbtab2empty_columns[i]
            new_tab = [sbtab[0].split('\n')[0]]
            #new_tab.append(sbtab[0].split('\n')[1])
            tab = sbtab[0].split('\n')[1:]
            for row in tab:
                new_row = []
                for j,element in enumerate(row.split('\t')):
                    if not j in eliminate_columns:
                        new_row.append(element)
                new_tab.append('\t'.join(new_row))
            new_tabs.append(('\n'.join(new_tab),sbtab[1]))

        return new_tabs
                

    def compartmentSBtab(self):
        '''
        Builds a Compartment SBtab.
        '''
        compartment  = [['!!SBtab SBtabVersion="1.0" Document="'+self.filename.rstrip('.xml')+'" TableType="Compartment" TableName="Compartment"'],['']]
        header       = ['!Compartment','!Name','!Size','!Unit','!SBOTerm']
        identifiers  = []
        column2ident = {}

        for comp in self.model.getListOfCompartments():
            value_row = ['']*len(header)
            value_row[0] = comp.getId()
            try: value_row[1] = comp.getName()
            except: pass
            try: value_row[2] = str(comp.getSize())
            except: pass
            try: value_row[3] = str(species.getUnits())
            except: pass
            if str(comp.getSBOTerm()) != '-1': value_row[5] ='SBO:%.7d'%comp.getSBOTerm()
            try:
                annot_tuples = self.getAnnotations(comp)
                for i,annotation in enumerate(annot_tuples):
                    if not ("!Identifiers:"+annotation[1]) in header:
                        identifiers.append(annotation[1])
                        column2ident[annotation[1]] = int(len(header)+1)
                        header.append('!Identifiers:'+annotation[1])
                        value_row.append('')
                    if annotation[1] in identifiers:
                        value_row[column2ident[annotation[1]]-1] = annotation[0]
            except: pass
            compartment.append(value_row)

        compartment[1] = header
        compartment_SB = compartment[0]
        for row in compartment[1:]:
            compartment_SB.append('\t'.join(row))
        compartment_SBtab = '\n'.join(compartment_SB)

        return [compartment_SBtab,'compartment']
        
    def compoundSBtab(self):
        '''
        Builds a Compound SBtab.
        '''
        compound = [['!!SBtab SBtabVersion="1.0" Document="'+self.filename.rstrip('.xml')+'" TableType="Compound" TableName="Compound"'],['']]
        header   = ['!Compound','!Name','!Location','!Charge','!IsConstant','!SBOTerm','!InitialConcentration','!hasOnlySubstanceUnits']
        identifiers  = []
        column2ident = {}

        for species in self.model.getListOfSpecies():
            value_row = ['']*len(header)
            value_row[0] = species.getId()
            value_row[1] = species.getName()
            try: value_row[2] = species.getCompartment()
            except: pass
            try: value_row[3] = str(species.getCharge())
            except: pass
            try: value_row[4] = str(species.getConstant())
            except: pass
            if str(species.getSBOTerm()) != '-1': value_row[5] ='SBO:%.7d'%species.getSBOTerm()
            try: value_row[6] = str(species.getInitialConcentration())
            except: pass
            try: value_row[7] = str(species.getHasOnlySubstanceUnits())
            except: pass
            try:
                annot_tuples = self.getAnnotations(species)
                for i,annotation in enumerate(annot_tuples):
                    if not ("!Identifiers:"+annotation[1]) in header:
                        identifiers.append(annotation[1])
                        column2ident[annotation[1]] = int(len(header)+1)
                        header.append('!Identifiers:'+annotation[1])
                        value_row.append('')
                    if annotation[1] in identifiers:
                        value_row[column2ident[annotation[1]]-1] = annotation[0]
            except: pass
            compound.append(value_row)

        compound[1] = header
        compound_SB = compound[0]
        for row in compound[1:]:
            compound_SB.append('\t'.join(row))
        compound_SBtab = '\n'.join(compound_SB)
            
        return [compound_SBtab,'compound']

    def eventSBtab(self):
        '''
        Builds an Event SBtab.
        '''
        if len(self.model.getListOfEvents()) == 0:
            return False
            
        event    = [['!!SBtab SBtabVersion="1.0" Document="'+self.filename.rstrip('.xml')+'" TableType="Event" TableName="Event"'],['']]
        header   = ['!Event','!Name','!Assignments','!Trigger','!SBOterm','!Delay','!UseValuesFromTriggerTime']
        identifiers  = []
        column2ident = {}

        for eve in self.model.getListOfEvents():
            value_row = ['']*len(header)
            value_row[0] = eve.getId()
            value_row[1] = eve.getName()
            if eve.getNumEventAssignments() > 1:
                try:
                    eas = eve.getListOfEventAssignments()
                    ea_entry = ''
                    for ea in eas:
                        var       = ea.getVariable()
                        ea_e      = ea.getMath()
                        ea_entry += var+' = '+libsbml.formulaToL3String(ea_e)+' | '
                    value_row[2] = ea_entry[:-3]
                except: pass
            else:
                try:
                    eas = eve.getListOfEventAssignments()
                    for ea in eas:
                        var          = ea.getVariable()
                        vr           = ea.getMath()
                        value_row[2] = var+' = '+libsbml.formulaToL3String(vr)
                except: pass
            try:
                trigger      = eve.getTrigger().getMath()
                value_row[3] = libsbml.formulaToL3String(trigger)
            except: pass
            if str(eve.getSBOTerm()) != '-1': value_row[4] ='SBO:%.7d'%eve.getSBOTerm()
            try: value_row[5] = str(eve.getDelay())
            except: pass
            try: value_row[6] = str(eve.getUseValuesFromTriggerTime())
            except: pass
            try:
                annot_tuples = self.getAnnotations(eve)
                for i,annotation in enumerate(annot_tuples):
                    if not ("!Identifiers:"+annotation[1]) in header:
                        identifiers.append(annotation[1])
                        column2ident[annotation[1]] = int(len(header)+1)
                        header.append('!Identifiers:'+annotation[1])
                        value_row.append('')
                    if annotation[1] in identifiers:
                        value_row[column2ident[annotation[1]]-1] = annotation[0]
            except: pass
            event.append(value_row)

        event[1] = header
        event_SB = event[0]
        for row in event[1:]:
            event_SB.append('\t'.join(row))
        event_SBtab = '\n'.join(event_SB)
            
        return [event_SBtab,'event']

    def ruleSBtab(self):
        '''
        Builds a Rule SBtab.
        '''
        if len(self.model.getListOfRules()) == 0:
            return False
            
        rule     = [['!!SBtab SBtabVersion="1.0" Document="'+self.filename.rstrip('.xml')+'" TableType="Rule" TableName="Rule"'],['']]
        header   = ['!Rule','!Name','!Formula','!Unit']
        identifiers  = []
        column2ident = {}

        for ar in self.model.getListOfRules():
            value_row = ['']*len(header)
            value_row[0] = ar.getId()
            value_row[1] = ar.getElementName()
            try:
                var          = ar.getVariable()
                vr           = ar.getMath()
                value_row[2] = var+' = '+libsbml.formulaToL3String(vr)
            except: pass
            try: value_row[3] = ar.getUnits()
            except: pass            
            try:
                annot_tuples = self.getAnnotations(ar)
                for i,annotation in enumerate(annot_tuples):
                    if not ("!Identifiers:"+annotation[1]) in header:
                        identifiers.append(annotation[1])
                        column2ident[annotation[1]] = int(len(header)+1)
                        header.append('!Identifiers:'+annotation[1])
                        value_row.append('')
                    if annotation[1] in identifiers:
                        value_row[column2ident[annotation[1]]-1] = annotation[0]
            except: pass
            rule.append(value_row)

        rule[1] = header
        rule_SB = rule[0]
        for row in rule[1:]:
            rule_SB.append('\t'.join(row))
        rule_SBtab = '\n'.join(rule_SB)
            
        return [rule_SBtab,'rule']

    def getAnnotations(self,element):
        '''
        Tries to extract an annotation from an SBML element.
        '''
        cvterms      = element.getCVTerms()
        annotation   = False
        urn          = False
        annot_tuples = []

        pattern2urn = {"CHEBI:\d+$":"obo.chebi",\
                       "C\d+$":"kegg.compound",\
                       "GO:\d{7}$":"obo.go",\
                       "((S\d+$)|(Y[A-Z]{2}\d{3}[a-zA-Z](\-[A-Z])?))$":"sgd",\
                       "SBO:\d{7}$":"biomodels.sbo",\
                       "\d+\.-\.-\.-|\d+\.\d+\.-\.-|\d+\.\d+\.\d+\.-|\d+\.\d+\.\d+\.(n)?\d+$":"ec-code",\
                       "K\d+$":"kegg.orthology",\
                       "([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\.\d+)?$":"uniprot"}
        for i in range(element.getNumCVTerms()):
            cvterm   = cvterms.get(i)
            for j in range(cvterm.getNumResources()):
                resource = cvterm.getResourceURI(j)
                for pattern in pattern2urn.keys():
                    try:
                        search_annot = re.search(pattern,resource)
                        annotation   = search_annot.group(0)
                        urn          = pattern2urn[pattern]
                        annot_tuples.append([annotation,urn])
                        break
                    except: pass
                    try:
                        search_annot = re.search('identifiers.org/(.*)/(.*)',resource)
                        annotation   = search_annot.group(2)
                        urn          = search_annot.group(1)
                        annot_tuples.append([annotation,urn])
                    except: pass

        return annot_tuples

    def reactionSBtab(self):
        '''
        Builds a Reaction SBtab.
        '''
        reaction     = [['!!SBtab SBtabVersion="1.0" Document="'+self.filename.rstrip('.xml')+'" TableType="Reaction" TableName="Reaction"'],['']]
        header       = ['!Reaction','!Name','!ReactionFormula','!Location','!Regulator','!KineticLaw','!SBOTerm','!IsReversible']
        identifiers  = []
        column2ident = {}

        for react in self.model.getListOfReactions():
            value_row = ['']*len(header)
            value_row[0] = react.getId()
            value_row[1] = react.getName()
            value_row[2] = self.makeSumFormula(react)
            try: value_row[3] = str(react.getCompartment())
            except: pass
            modifiers = react.getListOfModifiers()
            if len(modifiers)>1:
                modifier_list = ''
                for i,modifier in enumerate(modifiers):
                    if i != len(react.getListOfModifiers())-1: modifier_list += modifier.getSpecies() + '|'
                    else: modifier_list += modifier.getSpecies()
                value_row[4] = modifier_list
            elif len(modifiers)==1:
                for modifier in modifiers: value_row[4] = modifier.getSpecies()
            else: pass
            try:
                fm = react.getKineticLaw().getFormula()
                value_row[5] = fm.replace('\n','')
            except: pass
            if str(react.getSBOTerm()) != '-1': value_row[6] ='SBO:%.7d'%react.getSBOTerm()
            try: value_row[7] = str(react.getReversible())
            except: pass            
            try:
                annot_tuples = self.getAnnotations(react)
                for i,annotation in enumerate(annot_tuples):
                    if not ("!Identifiers:"+annotation[1]) in header:
                        identifiers.append(annotation[1])
                        column2ident[annotation[1]] = int(len(header))
                        header.append('!Identifiers:'+annotation[1])
                        value_row.append('')
                    if annotation[1] in identifiers:
                        value_row[column2ident[annotation[1]]] = annotation[0]
            except: pass
            reaction.append(value_row)

        reaction[1] = header
        reaction_SB = reaction[0]
        for row in reaction[1:]:
            reaction_SB.append('\t'.join(row))
        reaction_SBtab = '\n'.join(reaction_SB)
        
        return [reaction_SBtab,'reaction']

    def quantitySBtab(self):
        '''
        Builds a Quantity SBtab.
        '''
        pars = True
        if len(self.model.getListOfParameters()) == 0:
            pars = False
            for reaction in self.model.getListOfReactions():
                kinetic_law = reaction.getKineticLaw()
                if len(kinetic_law.getListOfParameters()) != 0: pars = True

        if not pars: return False
        
        quantity_SBtab = '!!SBtab SBtabVersion="1.0" Document="'+self.filename.rstrip('.xml')+'" TableType="Quantity" TableName="Quantity"\n!Quantity\t!Parameter:SBML:parameter:id\t!Value\t!Unit\t!Type'
        identifiers = []
        the_rows    = ''

        for reaction in self.model.getListOfReactions():
            kinetic_law = reaction.getKineticLaw()
            if kinetic_law:
                value_row   = ''
                for parameter in kinetic_law.getListOfParameters():
                    value_row += parameter.getId()+'_'+reaction.getId()+'\t'
                    value_row += parameter.getId()+'\t'
                    value_row += str(parameter.getValue())+'\t'
                    try: value_row += parameter.getUnits()+'\t'
                    except: value_row += '\t'
                    value_row += 'local parameter\t'
                    #quantity_SBtab += value_row
                    try:
                        (annotation,urn) = self.getAnnotation(parameter)
                        if urn in identifiers: value_row += annotation+'\n'
                        else: value_row += '\n'
                        if not "!Identifiers:" in quantity_SBtab:
                            quantity_SBtab += '\t'+'!Identifiers:'+urn
                            identifiers.append(urn)
                    except:
                        value_row += '\n'
                the_rows += value_row

        for parameter in self.model.getListOfParameters():
            value_row = parameter.getId()+'\t'
            value_row += parameter.getId()+'\t'
            value_row += str(parameter.getValue())+'\t'            
            try: value_row += parameter.getUnits()+'\t'
            except: value_row += '\t'            
            value_row += 'global parameter\t'
            try:
                (annotation,urn) = self.getAnnotation(parameter)
                if urn in identifiers: value_row += annotation+'\n'
                else: value_row += '\n'
                if not "!Identifiers:" in quantity_SBtab:
                    quantity_SBtab += '\t'+'!Identifiers:'+urn
                    identifiers.append(urn)
            except:
                value_row += '\n'
            the_rows += value_row

        quantity_SBtab += '\n'
        quantity_SBtab += the_rows
            
        return [quantity_SBtab,'quantity']

    def makeSumFormula(self,reaction):
        '''
        Generates the reaction formula of a reaction from the list of products and list of reactants.

        Parameters
        ----------
        reaction : libsbml object reaction
           Single reaction object from the SBML file.
        '''
        sumformula = ''
        id2name    = {}
        
        for species in self.model.getListOfSpecies():
            id2name[species.getId()] = species.getName()

        for i,reactant in enumerate(reaction.getListOfReactants()):
            if i != len(reaction.getListOfReactants())-1:
                if reactant.getStoichiometry() != 1.0:
                    sumformula += str(float(reactant.getStoichiometry())) + ' ' + reactant.getSpecies()+' + '
                else:
                    sumformula += reactant.getSpecies()+' + '
            else:
                if numpy.isnan(reactant.getStoichiometry()):
                    sumformula += '1 ' + reactant.getSpecies() + ' <=> '
                elif reactant.getStoichiometry() != 1.0:
                    sumformula += str(float(reactant.getStoichiometry())) + ' ' + reactant.getSpecies()+' <=> '
                else:
                    sumformula += reactant.getSpecies()+' <=> '
        
        if sumformula == '': sumformula += '<=> '
        
        for i,product in enumerate(reaction.getListOfProducts()):
            if i != len(reaction.getListOfProducts())-1:
                if product.getStoichiometry() != 1.0:
                    sumformula += str(float(product.getStoichiometry())) + ' ' + product.getSpecies()+' + '
                else:
                    sumformula += product.getSpecies()+' + '
            else:
                if numpy.isnan(product.getStoichiometry()):
                    sumformula += '1 ' + product.getSpecies() + ' <=> '
                elif product.getStoichiometry() != 1.0:
                    sumformula += str(float(product.getStoichiometry())) + ' ' + product.getSpecies()
                else:
                    sumformula += product.getSpecies()
            
        #if there is no product in the reaction (e.g. influxes), don't forget the tab
        #if len(reaction.getListOfProducts()) < 1:
        #    sumformula += '\t'

        return sumformula
        

if __name__ == '__main__':

    try: sys.argv[1]
    except:
        print('You have not provided input arguments. Please start the script by also providing an SBML file and an optional SBtab output filename: >python sbml2sbtab.py SBMLfile.xml Output')
        sys.exit()

    file_name  = sys.argv[1]
    try: output_name = sys.argv[2]+'.csv'
    except: output_name = file_name[:-4]+'.csv'

    reader     = libsbml.SBMLReader()
    sbml       = reader.readSBML(file_name)
    model      = sbml.getModel()
    Sbml_class = SBMLDocument(model,file_name)

    (sbtabs,warnings) = Sbml_class.makeSBtabs()

    #print warnings

    for sbtab in sbtabs:
        sbtab_name = output_name[:-4]+'_'+sbtab[1]+output_name[-4:]
        sbtab_file = open(sbtab_name,'w')
        sbtab_file.write(sbtab[0])
        sbtab_file.close()

    print('The SBtab file/s have been successfully written to your working directory or chosen output path.')
