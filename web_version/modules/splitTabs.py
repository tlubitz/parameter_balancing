#!/usr/bin/env python
import re, libsbml
import SBtab
import tablibIO
import misc

def checkTabs(document,filename,separator=None):
    '''
    this function checks, how many SBtab files are given by the user and saves it/them
    in a list, moreover stores the SBtab types in a dict linking to the SBtabs
    '''
    sbtabs     = []
    types      = []
    docs       = []
    tnames     = []
    type2sbtab = {}

    #if there are more than one SBtabs given in single files that might be comprised of several SBtabs:
    if len(document)>1:
        for single_document in document:
            #check for several SBtabs in one document
            document_rows = single_document.split('\n')
            tabs_in_document = getAmountOfTables(document_rows)
            if tabs_in_document > 1:
                sbtabs = splitDocumentInTables(document_rows)
            else: sbtabs = [document_rows]
            #generate SBtab class instance for every SBtab
            for sbtab in sbtabs:
                new_tablib_obj = tablibIO.importSetNew(sbtab,filename,separator)
                single_tab = SBtab.SBtabTable(new_tablib_obj,filename)
                type2sbtab[single_tab.table_type] = single_tab
                types.append(single_tab.table_type)
                docs.append(single_tab.table_document)
                try: tnames.append(single_tab.table_name)
                except: tnames.append('')
    #elif there is only one document given, possibly consisting of several SBtabs
    else:
        #check for several SBtabs in one document
        document_rows    = document[0].split('\n')
        tabs_in_document = getAmountOfTables(document_rows)
        if tabs_in_document > 1: sbtabs = splitDocumentInTables(document_rows)
        else: sbtabs = [document_rows]
        #generate SBtab class instance for every SBtab
        for sbtab in sbtabs:
            as_sbtab = '\n'.join(sbtab)
            new_tablib_obj = tablibIO.importSetNew(as_sbtab,filename,separator)
            single_tab = SBtab.SBtabTable(new_tablib_obj,filename)
            type2sbtab[single_tab.table_type] = single_tab
            types.append(single_tab.table_type)
            docs.append(single_tab.table_document)
            try: tnames.append(single_tab.table_name)
            except: tnames.append('')
            
    return sbtabs,types,docs,tnames

def getAmountOfTables(document_rows):
    '''
    counts the SBtab tables that are present in the document
    '''
    counter = 0
    for row in document_rows:
        if row.startswith('!!'):
            counter += 1

    return counter

def splitDocumentInTables(document_rows):
    '''
    if the document contains more than one SBtab, this function splits the document
    into the single SBtabs
    '''
    delimiter = misc.getDelimiter('\n'.join(document_rows))
    
    single_sbtab = [document_rows[0]]
    sbtab_list   = []
    for row in document_rows[1:]:
        if not row.startswith('!!'):
            splitrow = row.split(delimiter)
            if not row == '' and not row[0] == '' and not splitrow[0] == '' and not splitrow[0] == '':
                single_sbtab.append(row)
        else:
            sbtab_list.append(single_sbtab)
            single_sbtab = [row]
    sbtab_list.append(single_sbtab)

    return sbtab_list
