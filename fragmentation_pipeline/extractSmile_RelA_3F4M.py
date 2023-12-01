#!usr/bin/python

import pandas as pd
import pickle
import numpy as np
import os 
import sys
import rdkit
from rdkit import Chem
from rdkit.Chem import Recap
from rdkit.Chem.BRICS import BRICSDecompose
from rdkit.Chem import rdFMCS
import parser
import re

def removeNesting(elements):
    '''
    Removes the nested aspect of the input list
    input-spec: nested List or a simple list (max can do is 2 nested)
    output-spec: a set representation of the elements in the list

    TODO: use itertools to flatten the list instead
    '''
    flat=[]
    theSet = set() 
    for ele in elements:
        if type(ele) == list and len(ele) != 0:
            flat.extend(ele)
        elif type(ele) == list and len(ele) == 0 :
            continue
        elif ele != '':
            theSet.add(ele)
    for ele in flat:
        theSet.add(ele)
    return theSet

def matchingLigands(interactingAtoms,subStructureMatches):
    '''
    Matches the atom's ID from the subStrucure dictionary to the interacting Atoms;
    It is the interacting atom's ID that we are interested in; and we want to identify the
    ligands that have a single ID from that list.

    Input-spec: the first argument is the list of lists that contains the IDs of the interacting Atoms
                the second argumnent is the SET (IMPORTANT: SETS) that contains the matches
    Output-spec: it returns a SET of Ligands
    '''
    interactingSet = interactingAtoms
    matched=set()
    for key in subStructureMatches:
        if len(interactingSet.intersection(subStructureMatches[key])):
            matched.add(key)
    return matched


def translate(original,dics):
    '''
    This functions will translate between the atoms count that starts from 1 to etc to the
    atom's ID itself. It uses it with the SMILESTOPDB in the xml file
    input_spec: the first argument IS the subStructureMatches's IDs (dictionaries)
                the second argument IS the SMILESTOPDB dictionary
    output_spec: it DOES NOT RETURN A VALUE: it modifies the dictionary in place.
    '''
    for key in original:
        new=[]
        new=map(lambda x: dics[str(x+1)],original[key][0])
        original[key]=set(new)


def getSubStructIDs(ligand,fragments):
    '''
    obtains the fragment's atom IDs by matching them back to the original ligand's atom IDs
    input: ligand in smile format
            fragments in a list
    output-spec: a dictionary containing a key value pair in the format of the below
                "fragments:listofAtom'sIDs"
    '''
    mol_ligand = Chem.MolFromSmiles(ligand)
    theResult = {}
    for smallFragment in fragments:
        fragment_smile=smallFragment
        #params is needed to address the wildcards which may mess up the results
        params = Chem.AdjustQueryParameters()
        params.makeBondsGeneric = True
        params.makeDummiesQueries = True
        params.makeAtomsGeneric = True
        mol_fragment = Chem.MolFromSmiles(smallFragment)
        mol_fragment_adjusted= Chem.AdjustQueryProperties(mol_fragment, params)
        match_idx = mol_ligand.GetSubstructMatches(mol_fragment_adjusted)
        theResult[smallFragment]=match_idx
    return theResult

def getFragments(smile):
    '''
    This function returns the Fragments that is contained within the smile with the ReCap Decompose funcitonality
    input_spec: smiles presentation of a ligand
    output: a list of fragments
    '''
    m=Chem.MolFromSmiles(smile)
    expr = re.compile(r'[0-9]+\*')
    try:
        res = list(BRICSDecompose(m))
        newRes=[]
        oldRes={}
        for k in res:
            a = expr.sub('*',k)
            newRes.append(a)
            oldRes[a]=k
        print('the result of fragmentation: ',newRes)
        return (sorted(newRes),oldRes)
    except:
        raise ValueError

def getDictGenerated(string):
    '''
    it obtains the SMILESTOPDB_Conversion as a stirng and then stores it in the dictionary
    input-specs: the string with "key:value key2:value2" etc format
    output-specs: the dictionary that contains such representations.
    '''
    target = string.split(",")
    smileID=[]
    pdb_Id=[]
    for pairs in target:
        y = pairs.split(":")
        smileID.append(y[0])
        pdb_Id.append(y[1])
    zipObj = zip(smileID,pdb_Id)
    dictOfResults=dict(zipObj)
    return dictOfResults

def getCount(masterdict):
    countOfFragments={}
    iterator = 0
    for keys in masterdict:
        if float(masterdict[keys]['VINA']) < -6.0:
            iterator += 1
            fragments = masterdict[keys]["interactingLigands"]
            ligAndProDict=masterdict[keys]['LigToProDict']
            for fragment in fragments:
                if fragment in countOfFragments:
                    countOfFragments[fragment]['count'] += 1
                    countOfFragments[fragment]['ligand'].add(keys) 
                    for atom in masterdict[keys]['SubStructureMatches'][fragment]:
                        if atom in ligAndProDict:
                            theSet = countOfFragments[fragment]['protein']
                            countOfFragments[fragment]['protein']=theSet.union(ligAndProDict[atom])
                else:
                    information={}
                    information['count'] = 1
                    information['ligand'] = set()
                    information['ligand'].add(keys)
                    information['protein'] = set()
                    for atom in masterdict[keys]['SubStructureMatches'][fragment]:
                        if atom in ligAndProDict:
                            theSet = information['protein']
                            information['protein'] = theSet.union(ligAndProDict[atom])
                    countOfFragments[fragment]=information
        else:
            pass
    return countOfFragments 
def getAffinity(path,masterdict):
    '''
    This is real awkward; I just want to get the affinity but this FORCED ME to include a path to directory of ligands
    input-specs: path is to the ligand directories
                masterdict is the highest order of dictionary
    output-specs: it does not return anything
    '''
    for filename in os.listdir(path):
        tmpString=filename.split(".")
        name=tmpString[0]
        if len(name) <= 1:
            continue
        with open(path+filename) as file:
            for line in file:
                #Note: why is this here? is this still testing?
                #' '.join(line.split())
                #Note: use find command
                if line.find('VINA') != -1:
                    try:
                        masterdict[name]['VINA']=line.split('REMARK VINA RESULT:')[1].split()[0]
                    except:
                        print("Error: Could not get VINA score.")
def getFragToDom(fragToResDict,resToDomDict):
    theDict={}
    for frag in fragToResDict:
        domains=set()
        for protein in fragToResDict[frag]:
            print("{} vs {}".format(protein,resToDomDict.keys()))
            if protein in resToDomDict:
                domains.add(resToDomDict[protein])
                print("FOUND")
            else:
                domains.add("Unidentified Domain")
        theDict[frag]=domains

    return theDict

def getFragToRes(ligToProDict,ResToProDict,subStructureMatches,interactingLigands):
    '''
    To get the mapping from fragments to the sub structures
    input;  the Ligand Atom to Protein Atom Dictionary

            the Residue to Protein Dictionary

            the SubSutrcture Matches dictionary
    output: the dictionary mapping fragments to the Residues
    '''
    theDict={}
    for fragment in interactingLigands:
        theSet=subStructureMatches[fragment]
        proSet=set()
        for ele in theSet:
            if ele in ligToProDict:
                proSet=proSet.union(ligToProDict[ele])
        ResSet=set()

        for ele2 in proSet:
            if ele2 in ResToProDict:
                ResSet=ResSet.union(ResToProDict[ele2])
        theDict[fragment] = ResSet
    return theDict
            

def getAllInfoIntoDict(xmlDir,ligandDir):
    '''
    this method name is outdated; but it is best representation of what it dos
    input name: paths to both the xmlDir and ligandDir in that order
    output; a master dict file that would handle the info processing
    '''
    masterdict = {}
    iterator=0
    #resToDomDict=pickle.load( open('resToDom.p','rb'))
    for filename in os.listdir(xmlDir):
        details = {}
        tmpString = filename.split('-')
        print(tmpString[2])
        errors=[]
        try:
            dictsOfSmileToPDB = getDictGenerated(parser.getSMILE_PDBCONVERSION(xmlDir + filename))
            result =  parser.parse_ligand_indices(xmlDir + filename)
            interactingAtoms = result[0]      
            LigToProDict = result[1]
            residueDict=result[2]
            SMILE = parser.getSMILE(xmlDir + filename)
            Fragments = getFragments(SMILE)[0]
            originalFragments=getFragments(SMILE)[1]
            subStructureMatches = getSubStructIDs(SMILE, Fragments)
            translate(subStructureMatches, dictsOfSmileToPDB)
            matched = matchingLigands(interactingAtoms, subStructureMatches)
            fragToResDict=getFragToRes(LigToProDict,residueDict,subStructureMatches,matched) 
            #fragToDomDict=getFragToDom(fragToResDict,resToDomDict)
            details["SmiletoPDB"] = dictsOfSmileToPDB
            details["SMILE"] = SMILE
            details["Fragments"] = Fragments
            details["SubStructureMatches"] = subStructureMatches
            details["interactingLigands"] = matched
            details['LigToProDict']=LigToProDict
            details['ResidueToProtein']=result[2]
            details['FragToResDict'] = fragToResDict
            details['originalFragments']=originalFragments
        except Exception as e:
            print("ERROR: ",e)
            continue 
        #details['FragToDomDict']=fragToDomDict
        masterdict[tmpString[2]] = details
    getAffinity(ligandDir, masterdict)
#    addResidue(masterdict)
    return masterdict

def addResidue(masterdict):
    residueDict={}
    try:
        proDF=pd.read_csv('mapProToRes.csv')
    except:
        raise Exception("File not found")
    theDF=proDF.transpose()
    preDict=theDF.to_dict()
    for key in preDict:
        residueDict[preDict[key]['protein']]=preDict[key]['residueNumber']
    for ligand in masterdict:
        fragmentResidue={}
        for fragment in masterdict[ligand]['interactingLigands']:
            theSet=masterdict[ligand]['SubStructureMatches'][fragment]
            theLigToProDict=masterdict[ligand]['LigToProDict']
            proteinSet=set()
            for atom in theSet:
                if atom in theLigToProDict:
                    potentialSet=theLigToProDict[atom]
                    currentSet=proteinSet
                    proteinSet=currentSet.union(potentialSet)
            residueSet=set()
            for atom in proteinSet:
                residueSet.add(residueDict[int(atom)])
            
            fragmentResidue[fragment]=residueSet 
        masterdict[ligand]['FragmentResidues']=fragmentResidue
if __name__ == "__main__":
   
    if len(sys.argv) != 4:
        print('Needs 2 arguments. The directory of xml files, along with directory of ligands')
    else:
        ligandDir=sys.argv[2]
        xmlDir=sys.argv[1]
        tmpDir=sys.argv[3]
        mainPath=os.path.join(tmpDir,'main.csv')
        countPath=os.path.join(tmpDir,'output.csv')
        mainDir=getAllInfoIntoDict(xmlDir,ligandDir)
        picklePath=os.path.join(tmpDir,'main.pickle')
        with open(picklePath,'wb') as f:
            pickle.dump(mainDir,f)
        df4=pd.DataFrame(mainDir)
        df10=df4.reset_index()
        df11= df10.transpose()
        df11.rename({'Unnamed: 0':'Ligands'}, inplace=True ) 
        df11.rename({"index":"Ligands"},inplace=True)
        #df4=pd.DataFrame
        df11.to_csv(mainPath)
        
        #countDir=getCount(mainDir)
        #df5=pd.DataFrame(countDir)
        #df6=df5.transpose()
        #df6.to_csv(countPath)
 
    '''
        def2 = pd.DataFrame(masterdict)
        print(def2) 
        #testing out, not in production
        def2.to_excel('output.xlsx')

        df3= pd.read_excel('output.xlsx')
        print(df3)
        newdict=df3.to_dict()
        print(newdict)
    '''
