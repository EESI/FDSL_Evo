from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import QED
from rdkit.Chem import Descriptors
from identifyAtomOnLigand import getBestCarbon
import random
import argparse
import os
import sys
import re
import pandas as pd
import pickle

def replaceHydrogenWithDummy(mol, row):
    indexToReplace = None
    for bond in mol.GetAtomWithIdx(row).GetBonds():
        if (bond.GetEndAtom().GetAtomicNum() == 1):
            indexToReplace = bond.GetEndAtom().GetIdx()
        
    if indexToReplace != None:
        mol.GetAtomWithIdx(indexToReplace).SetAtomicNum(0)

    return mol

def select_random_carbon_with_hydrogen(mol):
    carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            # Check if the carbon has a hydrogen attached
            if any(nb.GetAtomicNum() == 1 for nb in atom.GetNeighbors()):
                carbons.append(atom.GetIdx() + 1)

    if not carbons:
        raise ValueError("No carbon atoms with a hydrogen attached found in the molecule.")

    random_carbon_row = random.choice(carbons)
    return random_carbon_row

def findDummy(MOL):
    """
    finds the dummy atom "*" and the proper attachment point, cname.
    """
    conn_atoms = [a.GetIdx() for a in MOL.GetAtoms() if a.GetAtomicNum() == 0] #get the index for all atomic number == 0 "atoms"
    neighbors = [MOL.GetAtomWithIdx(x).GetNeighbors()[0].GetIdx() for x in conn_atoms] #find neighboring atom (where the connection will actually be made)
    print("NEIGHBORS",neighbors[0])
    return neighbors[0]


def removeDummies(MOL):
    try:
        listDummies = []
        print("BEFORE REMOVING DUMMIES")
        print(Chem.MolToMolBlock(MOL))
        for a in MOL.GetAtoms():
            if a.GetSymbol() == '*' or a.GetSymbol() == 'R':
                listDummies.append(a.GetIdx())
        for i in range(len(listDummies)):
            MOL.RemoveAtom(listDummies[i]-i)
        #MOL.RemoveAtom([a.GetIdx() for a in MOL.GetAtoms() if a.GetAtomicNum() == 0][0])
        print("AFTER REMOVING DUMMIES")
        print(Chem.MolToMolBlock(MOL))
    except Exception as e: #if no dummies
        print("REMOVE DUMMIES ERROR", e)
    return MOL

def MMFF94opt(mol):
    Chem.SanitizeMol(mol) 
    mol=Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=739)
    AllChem.MMFFOptimizeMolecule(mol) #opt because spawned fragments are at origin, not proper 3D place on existing mol.
    final_mol = Chem.RemoveHs(mol) #not entirely effective. Need openbabel's help to remove all H with -d flag.
    return final_mol

def mergeFrags(mol, frag, molIndex):
    try:
        merged = Chem.RWMol(Chem.CombineMols(mol, frag))
        attachmentPoint = findDummy(merged)
        merged.AddBond(molIndex,attachmentPoint,Chem.rdchem.BondType.SINGLE)

        merged = removeDummies(merged)
    except Exception as e:
        print("MERGE FRAGS ERROR", e)
    try:
        opt = MMFF94opt(merged)
    except Exception as e:
        print("MMFF94 ERROR", e)
    return opt


def extract_atom_range(mol_file, start, stop, output_file):
    print("TYPE:", type(mol_file))
    try:
        mol = Chem.MolFromMolFile(str(mol_file))
    except Exception as e:
        print("ERROR", e)

    if mol is None:
        print("Error reading MOL file:", mol_file)
        return

    num_atoms = mol.GetNumAtoms()
    print("NUM ATOMS", num_atoms)


    if start < 0 or start >= num_atoms or stop < 0 or stop >= num_atoms or start >= stop:
        print("Invalid start and stop indices.")
        return

    atom_range = set(range(start, stop+1))
    extracted_mol = Chem.RWMol(mol)

    for idx in reversed(range(num_atoms)):
        if idx not in atom_range:
            extracted_mol.RemoveAtom(idx)

    extracted_mol.UpdatePropertyCache(strict=False)
    extracted_mol = extracted_mol.GetMol()
    return extracted_mol 


def get_random_mol_from_dict(dict_pickle_path, AAID):
    # Read the pickle file into a dict
    try:
        with open(dict_pickle_path, 'rb') as handle:
            AADict = pickle.load(handle)

        print("AAID",AAID)

        # Select a random fragment
        random_fragment = random.choice(list(AADict[str(AAID)]))

        print("BEFORE SUB",random_fragment)

        # Replace the "[#*]" pattern with "*"
        random_fragment = re.sub(r'\[\d+\*\]', '*', random_fragment)

        print("AFTER SUB", random_fragment)

        # Convert SMILES fragment to 3D mol object
        mol = Chem.MolFromSmiles(random_fragment)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)      
    return mol

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("-x", "--xml", help="XML from plip output")
    parser.add_argument("-p", "--pronatedPDB", help="Pronated Merged PDB file from merge-it.py including ligand and receptor")
    parser.add_argument("-m", "--mergedMol", help="Mol File from obabel. Cleans merged output pdb")
    parser.add_argument("-g", "--gen", help="Generation currently on")
    parser.add_argument("-r", "--rank", help="MPI Rank")
    parser.add_argument("-a","--aadictionary", help="Path to AADictionary.pickle outputted from generateFragmentBank.sh")


    args = parser.parse_args()

    attempts = 0
    while True and attempts < 10:    
        attempts+=1
        try:
            picklePath = args.aadictionary 
            protonatedPDBPath = args.pronatedPDB
            molPath = args.mergedMol
            gen = args.gen
            rank = args.rank

            AAID, bestCarbon, start, stop= getBestCarbon(protonatedPDBPath, xmlPath, attempts)

            originalLigandMol = extract_atom_range(molPath, start-1, stop-1, "Test.mol")

            print("PICKLE PATH:", picklePath)
            mergedMol = mergeFrags(originalLigandMol, get_random_mol_from_dict(picklePath, AAID[:-1]), bestCarbon - start)
            print("FINISHED MERGE")
            
            Chem.MolToPDBFile(mergedMol, os.getenv('TMPDIR') + "/pdbForPLIP/" + rank + "/fragment_addition_result_" + gen + "_" + rank + ".pdb")
            print("Rank: " + str(rank)+ " FINISHED ADDITION. STORED AT" + os.getenv('TMPDIR') + "/pdbForPLIP/" + str(rank) + "/fragment_addition_result_" + str(gen) + "_" + str(rank) + ".pdb")
            break
        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)            
            

    