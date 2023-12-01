from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import QED
from rdkit.Chem import Descriptors
import sys
import math
import xml.etree.ElementTree as ET


#gets average coordinates of target AA
def getAACoordinates(pdbPath, aaIdentity, aaPosition):
    pdbFile = open(pdbPath,"r")
    xPos, yPos, zPos, count = 0, 0, 0, 0
    for atom in pdbFile.readlines():
        splitString = atom.split()
        if splitString[0] == "ATOM" and splitString[3] == aaIdentity and splitString[5] == aaPosition:
            xPos += float(splitString[6])
            yPos += float(splitString[7])
            zPos += float(splitString[8])
            count+=1
    
    pdbFile.close()
    return xPos/count,yPos/count,zPos/count
    
#gets closest carbon on ligand to target amino acid in protein
def getClosestCarbon(pdbPath, aaIdentityAndPosition, aaCoordinates):
    pdbFile = open(pdbPath,"r")
    lowestDistance = sys.maxsize
    closestCarbon = None

    for atom in pdbFile.readlines():
        splitString = atom.split()
        if splitString[0] == "HETATM" and (splitString[11] == "C" or splitString[11] == "N"):
            distance = math.sqrt((float(splitString[6]) - aaCoordinates[0])**2 
                                + (float(splitString[7]) - aaCoordinates[1])**2 
                                + (float(splitString[8]) - aaCoordinates[2])**2)
            if distance < lowestDistance:
                closestCarbon = splitString
                lowestDistance = distance
    
    pdbFile.close()
    return closestCarbon, aaIdentityAndPosition, lowestDistance

#determines closet ligand carbon to each protein amino acid
def findClosestCarbonToEveryAA(pdbPath):
    pdbFile = open(pdbPath,"r")
    completedAAs = []
    closestCarbonsAndDistances = []
    start = None
    stop = None

    for atom in pdbFile.readlines():
        splitString = atom.split()
        if splitString[0] == "ATOM" and splitString[5] not in completedAAs:
            closestCarbonsAndDistances.append(getClosestCarbon(pdbPath, splitString[3]+splitString[5], getAACoordinates(pdbPath, splitString[3], splitString[5])))
            completedAAs.append(splitString[5])
        if splitString[0] == "HETATM" and splitString[11] != "H":
            if start == None:
                start = splitString[1]
            stop = splitString[1]
    pdbFile.close()
    closestCarbonsAndDistances.sort(key = lambda x: x[2])
    return closestCarbonsAndDistances, int(start), int(stop)

#returns carbon closest located closest to any amino acid in protein
def getBestCarbon(pdbPath, xmlPath, attempts):
    listClosestCToAA, start, stop = findClosestCarbonToEveryAA(pdbPath)
    xmlTree = ET.parse(xmlPath)
    xmlRoot = xmlTree.getroot()
    bs_residues = xmlRoot.findall('.//bs_residue')
    print(bs_residues)
    for possibleAtom in listClosestCToAA:
        for bs_residue in bs_residues:
            id_value = bs_residue.text
            if id_value[:-1] == possibleAtom[1][3:]:
                print("NEW MATCH")
                print(id_value, int(possibleAtom[0][1]), start, stop)
                attempts -= 1
                if attempts == 0:
                    return id_value, int(possibleAtom[0][1]), start, stop
                break
if __name__ == "__main__":
    pass