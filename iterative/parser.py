import pandas as pd 
import xmltodict 
import ast 
import os
from json import loads, dumps
from collections import OrderedDict
import re
from ast import literal_eval
from  itertools import chain 


def json_to_list(input_ordered_dict):
    return loads(dumps(input_ordered_dict))

def addValToDicOfSets(theDict,key,value):
    theDict.setdefault(key, set())
    if isinstance(value,set):
        for i in value:
            theDict[key].add(i)
    else:
        theDict[key].add(value)
def getDataFromMultipleDict(theList,atomSet,LigAndProDict):
    for i in theList:
        prot=set()
        ligs=set()
        if isinstance(i['prot_idx_list']['idx'],list):
            for j in i['prot_idx_list']['idx']:   
                prot.add(j['#text']) 
        else:
            prot.add(i['prot_idx_list']['idx']['#text'])

        if isinstance(i['lig_idx_list']['idx'],list):
            for j in i['lig_idx_list']['idx']:   
                atomSet.add(j['#text'])
                ligs.add(j['#text'])
        else:
            atomSet.add(i['lig_idx_list']['idx']['#text'])
            ligs.add(i['lig_idx_list']['idx']['#text'])
        for lig in ligs:
            addValToDicOfSets(LigAndProDict,lig,prot)


def getDataFromSingleDict(theList,atomSet,LigAndProDict):
    prot=set()
    ligs=set()
    if isinstance(theList['prot_idx_list']['idx'],list):
        for j in theList['prot_idx_list']['idx']:   
            prot.add(j['#text']) 
    else:
        prot.add(theList['prot_idx_list']['idx']['#text'])

    if isinstance(theList['lig_idx_list']['idx'],list):
        for j in theList['lig_idx_list']['idx']:   
            atomSet.add(j['#text'])
            ligs.add(j['#text'])
    else:
        ligs.add(theList['lig_idx_list']['idx']['#text'])
        atomSet.add(theList['lig_idx_list']['idx']['#text'])
    for lig in ligs:
        addValToDicOfSets(LigAndProDict,lig,prot)

def addValToResDict(residueName,protein,ProToResDict):
    '''
    if residueName in ProToResDict:
        ProToResDict[residueName].add(protein)
    else:
        protSet=set()
        protSet.add(protein)
        ProToResDict[residueName]=protSet 
    '''
    if protein in ProToResDict:
        ProToResDict[protein].add(residueName)
    else:
        resSet=set()
        resSet.add(residueName)
        ProToResDict[protein]=resSet
        
def getDataFromMultipleDictForRes(theList,ProToResDict):
    for i in theList:
        residueName=i['restype']+' '+i['resnr']
        if isinstance(i['prot_idx_list']['idx'],list):
            for j in i['prot_idx_list']['idx']:   
                addValToResDict(residueName,j['#text'],ProToResDict)
        else:
            addValToResDict(residueName,i['#text'],ProToResDict)
            prot.add(i['prot_idx_list']['idx']['#text'])

def getDataFromSingleDictForRes(theList,ProToResDict):
    residueName=theList['restype'] + ' ' + theList['resnr']
    if isinstance(theList['prot_idx_list']['idx'],list):
        for j in theList['prot_idx_list']['idx']:   
            addValToResDict(residueName,j['#text'],ProToResDict)
    else:
        addValToResDict(residueName,theList['prot_idx_list']['idx']['#text'])
def getSMILE_PDBCONVERSION(path):
    xml_string = open(path).read()
    smile_to_pdb_start = re.search("<smiles_to_pdb>", xml_string).start()
    smile_to_pdb_end = re.search("</smiles_to_pdb>", xml_string).start() + len("</smiles_to_pdb>")
    smile_to_pdb = xml_string[smile_to_pdb_start:smile_to_pdb_end]
    smile_to_pdb = json_to_list(xmltodict.parse(smile_to_pdb))["smiles_to_pdb"]
    return smile_to_pdb

def getSMILE(path):
    xml_string = open(path).read()
    smile_start = re.search("<smiles>",xml_string).start()
    smile_end = re.search("</smiles>",xml_string).start() + len("</smiles>")
    smile = xml_string[smile_start:smile_end]
    smile = json_to_list(xmltodict.parse(smile))["smiles"]
    return smile


def parse_ligand_indices(xml_path):
    xml_string = open(xml_path).read()
    smile_to_pdb_start = re.search("<smiles_to_pdb>", xml_string).start()
    smile_to_pdb_end = re.search("</smiles_to_pdb>", xml_string).start() + len("</smiles_to_pdb>")
    smile_to_pdb = xml_string[smile_to_pdb_start:smile_to_pdb_end]
    smile_to_pdb = json_to_list(xmltodict.parse(smile_to_pdb))["smiles_to_pdb"]
    smile_index_list = []
    atomSet = set()
    LigAndProDict={}
    protSet= set()
    ProToResDict={}
    for i in smile_to_pdb.split(","):
        smile_index_list.append(i.split(":")[1])
    try:
        hydrophobic_interactions_start = re.search("<hydrophobic_interactions>", xml_string).start()
        hydrophobic_interactions_end = re.search("</hydrophobic_interactions>", xml_string).start() + len("</hydrophobic_interactions>")
        hydrophobic_interactions = xml_string[hydrophobic_interactions_start:hydrophobic_interactions_end]
        hydrophobic_interactions = json_to_list(xmltodict.parse(hydrophobic_interactions))
        hydrophobic_interactions_list = []
        theDict=json_to_list(hydrophobic_interactions['hydrophobic_interactions']['hydrophobic_interaction'])
        if isinstance(theDict,list):
            for k in theDict: 
                residueName=k['restype'] + " "+ k['resnr']
                addValToResDict(residueName,k['protcarbonidx'],ProToResDict)
                atomSet.add(k["ligcarbonidx"])
                addValToDicOfSets(LigAndProDict,k['ligcarbonidx'],k['protcarbonidx'])
        else:
            residueName=theDict['restype'] + " "+ theDict['resnr']
            addValToResDict(residueName,theDict['protcarbonidx'],ProToResDict)
            atomSet.add(theDict["ligcarbonidx"])
            addValToDicOfSets(LigAndProDict,theDict['ligcarbonidx'],theDict['protcarbonidx'])

    except Exception as e:
        hydrophobic_interactions = ""
        hydrophobic_interactions_list = []
    try:
        hydrogen_bonds_start = re.search("<hydrogen_bonds>", xml_string).start()
        hydrogen_bonds_end = re.search("</hydrogen_bonds>", xml_string).start() + len("</hydrogen_bonds>")
        hydrogen_bonds = xml_string[hydrogen_bonds_start:hydrogen_bonds_end]
        hydrogen_bonds =  json_to_list(xmltodict.parse(hydrogen_bonds))["hydrogen_bonds"]["hydrogen_bond"]
        hydrogen_bonds_idx_list = []
        if isinstance(hydrogen_bonds,list):
            for i in hydrogen_bonds:
                residueName= (i['restype'] + ' '+ i['resnr'])
                if i["acceptoridx"] in smile_index_list:
                    atomSet.add(i["acceptoridx"])
                    addValToResDict(residueName,i['donoridx'],ProToResDict)
                    addValToDicOfSets(LigAndProDict,i['acceptoridx'],i['donoridx'])
                elif i["donoridx"] in smile_index_list:
                    addValToDicOfSets(LigAndProDict,i['donoridx'],i['acceptoridx'])
                    addValToResDict(residueName,i['acceptoridx'],ProToResDict)
                    atomSet.add(i["donoridx"])
        else:
            residueName = hydrogen_bonds['restype'] + ' ' + hydrogen_bonds['resnr']
            if hydrogen_bonds['acceptoridx'] in smile_index_list:
                addValToResDict(residueName,hydrogen_bonds['donoridx'],ProToResDict)
                atomSet.add(hydrogen_bonds["acceptoridx"])
                addValToDicOfSets(LigAndProDict,hydrogen_bonds['acceptoridx'],hydrogen_bonds['donoridx'])
            elif hydrogen_bonds['donoridx'] in smile_index_list:
                atomSet.add(hydrogen_bonds["donoridx"])
                addValToResDict(residueName,hydrogen_bonds['acceptoridx'],ProToResDict)
                addValToDicOfSets(LigAndProDict,hydrogen_bonds['donoridx'],hydrogen_bonds['acceptoridx'])

    except Exception as e:
        hydrogen_bonds = ""
        hydrogen_bonds_idx_list = []

    try:
        salt_bridges_start = re.search("<salt_bridges>", xml_string).start()
        salt_bridges_end = re.search("</salt_bridges>", xml_string).start() + len("</salt_bridges>")
        salt_bridges = xml_string[salt_bridges_start:salt_bridges_end]
        salt_bridges = json_to_list(xmltodict.parse(salt_bridges))
        salt_bridges_idx_list = []
        theList=json_to_list(salt_bridges['salt_bridges']['salt_bridge'])
        if isinstance(theList,dict):
            getDataFromSingleDict(theList,atomSet,LigAndProDict)    
            getDataFromSingleDictForRes(theList,ProToResDict)
        else:
            getDataFromMultipleDictForRes(theList,ProToResDict)
            getDataFromMultipleDict(theList,atomSet,LigAndProDict)
        '''
        for i in theList:
            print(i)
            prot=set()
            ligs=set()
        '''
    except:
        salt_bridges = ""
        salt_bridges_idx_list = []
    try:
        pi_stacks_start = re.search("<pi_stacks>", xml_string).start()
        pi_stacks_end = re.search("</pi_stacks>", xml_string).start() + len("</pi_stacks>")
        pi_stacks = xml_string[pi_stacks_start:pi_stacks_end]
        pi_stacks = json_to_list(xmltodict.parse(pi_stacks))
        theList= json_to_list(pi_stacks['pi_stacks']['pi_stack'])
        if isinstance(theList,dict):
            getDataFromSingleDict(theList,atomSet,LigAndProDict)    
            getDataFromSingleDictForRes(theList,ProToResDict)
        else:
            getDataFromMultipleDictForRes(theList,ProToResDict)
            getDataFromMultipleDict(theList,atomSet,LigAndProDict)

    except Exception as e:
        pi_stacks = ""
        pi_stacks_idx_list = []
    try:
        pi_cation_interactions_start = re.search("<pi_cation_interactions>", xml_string).start()
        pi_cation_interactions_end = re.search("</pi_cation_interactions>", xml_string).start() + len("</pi_cation_interactions>")
        pi_cation_interactions = xml_string[pi_cation_interactions_start:pi_cation_interactions_end]
        pi_cation_interactions =  json_to_list(xmltodict.parse(pi_cation_interactions))
        pi_cation_idx_list = []
        theList=json_to_list(pi_cation_interactions['pi_cation_interactions']['pi_cation_interaction'])
        if isinstance(theList,dict):
            getDataFromSingleDict(theList,atomSet,LigAndProDict)    
            getDataFromSingleDictForRes(theList,ProToResDict)
        else:
            getDataFromMultipleDictForRes(theList,ProToResDict)
            getDataFromMultipleDict(theList,atomSet,LigAndProDict)
        '''
            for i in theList:
                prot=[]
                ligs=[]
                if isinstance(i['prot_idx_list']['idx'],list):
                    for j in i['prot_idx_list']['idx']:   
                        prot.append(j['#text']) 
                else:
                    prot.append(i['prot_idx_list']['idx']['#text'])

                if isinstance(i['lig_idx_list']['idx'],list):
                    for j in i['lig_idx_list']['idx']:   
                        atomSet.add(j['#text'])
                        ligs.append(j['#text'])
                else:
                    atomSet.add(i['lig_idx_list']['idx']['#text'])
                    ligs.append(i['lig_idx_list']['idx']['#text'])
                for lig in ligs:
                    LigAndProDict[lig]=prot

            lig_idx_list = []
            for j in i["lig_idx_list"]["idx"]: 
                lig_idx_list.append(j["#text"])
                atomSet.add(j["#text"])
            for k in i['prot_idx_list']['idx']:
                pro_idx_list.append(k['#text'])
            for ele in i:
                LigAndProDict[ele]=pro_idx_list
            pi_cation_idx_list.append(lig_idx_list)
        '''
    except Exception as e:
        pi_cation_interactions = ""
        pi_cation_idx_list = []
    try:
        halogen_bonds_start = re.search("<halogen_bonds>", xml_string).start()
        halogen_bonds_end = re.search("</halogen_bonds>", xml_string).start() + len("</halogen_bonds>")
        halogen_bonds = xml_string[halogen_bonds_start:halogen_bonds_end]
        halogen_bonds = json_to_list(xmltodict.parse(halogen_bonds))    
    #     Is donor or acceptor of interest
        a = json_to_list(halogen_bonds['halogen_bonds']['halogen_bond']) 
        if isinstance(a,dict):
            residueName = a['restype'] + ' ' + a['resnr']
            if a["don_idx"] in smile_index_list:
                addValToResDict(residueName,a['acc_idx'],ProToResDict)
                atomSet.add(a["don_idx"])
                addValToDicOfSets(LigAndProDict,a['don_idx'],a['acc_idx'])
            elif a["acc_idx"] in smile_index_list:
                addValToResDict(residueName,a['don_idx'],ProToResDict)
                addValToDicOfSets(LigAndProDict,a['acc_idx'],a['don_idx'])
                atomSet.add(a["acc_idx"])
        else:
            for i in a:
                residueName = i['restype'] + ' ' + i['resnr']
                if i["don_idx"] in smile_index_list:
                    atomSet.add(i["don_idx"])
                    addValToResDict(residueName,i['acc_idx'],ProToResDict)

                    addValToDicOfSets(LigAndProDict,i['don_idx'],i['acc_idx'])
                elif i["acc_idx"] in smile_index_list:
                    addValToResDict(residueName,i['don_idx'],ProToResDict)
                    addValToDicOfSets(LigAndProDict,i['acc_idx'],i['don_idx'])
                    atomSet.add(i["acc_idx"])
 
    except Exception as e:
        halogen_bonds = ""
    try:
        metal_complexes_start = re.search("<metal_complexes>", xml_string).start()
        metal_complexes_end = re.search("</metal_complexes>", xml_string).start() + len("</metal_complexes>")
        metal_complexes = xml_string[metal_complexes_start:metal_complexes_end]
        metal_complexes = json_to_list(xmltodict.parse(metal_complexes))
    #     Is donor or acceptor of interest
        atomSet.add('metalic Complxes')
    except:
        metal_complexes = ""

#    return [hydrophobic_interactions_list, hydrogen_bonds_idx_list, salt_bridges_idx_list,  pi_stacks_idx_list,pi_cation_idx_list, halogen_bonds, metal_complexes  ]
    return (atomSet,LigAndProDict,ProToResDict)

