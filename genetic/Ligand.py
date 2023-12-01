from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import QED
from rdkit.Chem import Descriptors
import os

class Ligand:
    def __init__(self, bitList):
        self.bitList = bitList
        self.SMILES = None
        self.VINA = 0
        self.MW = None
        self.MOL = None
        self.QED = 0

    def getBitList(self):
        return self.bitList

    def __str__(self):
        return str(self.bitList)

    def getSMILES(self):
        return self.SMILES

    def setSMILES(self, value):
        self.SMILES = value

    def getVINA(self):
        return self.VINA
    
    def setVINA(self, value):
        self.VINA = value

    def getMW(self):
        return Descriptors.MolWt(self.getMOL())
   
    def setMW(self, value):
        self.MW = value

    def getMOL(self):
        return self.MOL
    
    def setMOL(self, mol):
        self.MOL = mol

    def getQED(self):
        return self.QED
    
    def setQED(self):
        self.QED = QED.qed(self.getMOL())

    #Manually add source ligand VINA/QED means and standard deviations if interested in multiobjective function
    def getMULTI(self):
        VINA_MEAN = -6.320103075178621
        VINA_STD = 0.9216458974454603
        QED_MEAN = 0.7110247681253329
        QED_STD = 0.12908877169218852
        return (self.getVINA() - VINA_MEAN)/VINA_STD - (self.getQED() - QED_MEAN)/QED_STD
    

    def __lt__(self, other):
        return self.VINA < other.getVINA()

    def __gt__(self, other):
        return self.VINA > other.getVINA()

    def __eq__(self, other):
        if self.getSMILES() == other.getSMILES() or self.getBitList() == other.getBitList():
            return True
        return False        

    def writeLigand(self, rank, mol2D, name):
        if int(rank) > 0:
            smiles = Chem.MolToSmiles(mol2D)
            self.setSMILES(smiles)
            #replaces dummy atoms with hydrogens (rarely failes)
            for a in mol2D.GetAtoms():
                if a.GetSymbol() == '*':
                    a.SetAtomicNum(1)
            try:
                molFromSmilesWithH = Chem.AddHs(mol2D)
                AllChem.EmbedMolecule(molFromSmilesWithH, randomSeed=739)   #add standardized seed
                AllChem.MMFFOptimizeMolecule(molFromSmilesWithH)
                embeddedMolecule = Chem.RemoveHs(molFromSmilesWithH)
                self.setMOL(embeddedMolecule)
                self.setQED()
                TMPDIRNAME = os.getenv('TMPDIR')
                Chem.MolToMolFile(embeddedMolecule, TMPDIRNAME+"/molFiles/"+name+".mol")
                return 1
            except Exception as e:
                print(e)
                return -1
    
    