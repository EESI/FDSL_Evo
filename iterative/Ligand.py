from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import QED
from rdkit.Chem import Descriptors
import os

class Ligand:
    def __init__(self, bitList, SMILES):
        self.bitList = bitList
        self.SMILES = SMILES
        self.VINA = 0
        self.MW = None
        self.MOL = None
        self.PREVLIG = None
        self.QED = 0
        self.atomSet = None
        self.ligToProDict = None
        self.proToResDict = None
        self.AAInteractions = None
        self.FragEnds = 0
        self.TASK = None
        self.GEN = 0
        self.RANK = None
        self.finished = False

    def getBitList(self):
        return self.bitList

    def __str__(self):
        return str(self.getSMILES())

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
   
    def setMW(self, df):
        numFragEnds = 0
        self.MW = Descriptors.MolWt(self.getMOL())
        self.FragEnds = numFragEnds

    def getMOL(self):
        return self.MOL
    
    def setMOL(self, mol):
        self.MOL = mol

    def getQED(self):
        return QED.qed(self.getMOL())
    
    def setQED(self):
        self.QED = QED.qed(self.getMOL())

    def getAtomSet(self):
        return self.atomSet
        
    def setAtomSet(self, value):
        self.atomSet = value

    def getLigToProDict(self):
        return self.ligToProDict

    def setLigToProDict(self, value):
        self.ligToProDict = value
    
    def getProToResDict(self):
        return self.proToResDict

    def setProToResDict(self, value):
        self.proToResDict = value

    def getAAInteractions(self):
        return self.AAInteractions

    def setAAInteractions(self, value):
        self.AAInteractions = value
    
    def getPreviousLigand(self):
        return self.PREVLIG

    def setPreviousLigand(self, ligand):
        self.PREVLIG = ligand

    def setTASK(self, value):
        self.TASK = value
    
    def getTASK(self):
        return self.TASK

    def getGEN(self):
        return self.GEN
    
    def setGEN(self, value):
        self.GEN = value

    def getRANK(self):
        return self.RANK

    def setRANK(self, value):
        self.RANK = value

    def setFinished(self, value):
        self.finished = value

    def isFinished(self):
        return self.finished

    #Manually add source ligand VINA/QED means and standard deviations if interested in multiobjective function
    def getMULTI(self):
        VINA_MEAN = -6.320103075178621
        VINA_STD = 0.9216458974454603
        QED_MEAN = 0.7110247681253329
        QED_STD = 0.12908877169218852
        return (self.getVINA() - VINA_MEAN)/VINA_STD - (self.getQED() - QED_MEAN)/QED_STD

    def getFragmentEnds(self):
        return self.FragEnds
    

    def __lt__(self, other):
        return self.VINA < other.getVINA()

    def __gt__(self, other):
        return self.VINA > other.getVINA()

    def __eq__(self, other):
        if other == None:
            if self.getSMILES() == None:
                return True
            return False
        if self.getSMILES() == other.getSMILES() or self.getBitList() == other.getBitList():
            return True
        return False        

    def reset(self, currentMol, df):
        self.setMOL(currentMol)
        self.setVINA(None)
        self.setSMILES(Chem.MolToSmiles(self.getMOL()))
        self.setMW(df)

    def writeLigand(self, rank, mol2D, name):
        if int(rank) > 0:
            if self.getMOL() == None:
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
            else:
                Chem.MolToMolFile(self.getMOL(), os.getenv('TMPDIR')+"/molFiles/"+name+".mol")

    def isEmpty(self):
        if self.getSMILES() == None:
            return True
        return False