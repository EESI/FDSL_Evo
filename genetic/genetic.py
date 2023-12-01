from numpy.random import randint
from numpy.random import rand
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import BRICS
from rdkit.Chem import AllChem
import random
import re
import subprocess
import os
from os import walk
from mpi4py import MPI
import matplotlib.pyplot as plt
from Ligand import Ligand
import copy
import argparse
import math


#Builds ligands using BRICS and tracks Build Fails   
def maxVina(ligand, df, name):
    if int(rank) > 0:
        random.seed(20)
        totalSum = 0
        smilesFragList = []
        for i in ligand.getBitList():
            binaryString = ''
            for j in i:
                binaryString+=str(j)
            index = int(binaryString, 2)
            smilesFrag = df.at[index, "Fragment"]
            smilesFragList.append(smilesFrag)
        molFrags = [Chem.MolFromSmiles(f) for f in smilesFragList]

        try:
            BRICSObject = BRICS.BRICSBuild(molFrags, maxDepth=1, scrambleReagents=False, onlyCompleteMols=True) #onlyCompleteMols=False to keep frags
            result = list(BRICSObject)[0]
        except:
            result = None
        if result != None:
            Chem.SanitizeMol(result)
            return ligand.writeLigand(rank, result, name)
        return 0  

#evalutes population through Autodock VINA
def evaluate(objective, pop, gen, df, ligandsTested, config_path, proteinPDB_path):
    failedBuild = 0
    failedRemoval = 0

    try:
        subprocess.run("mkdir $TMPDIR/molFiles > /dev/null 2>&1", shell = True)  
    except:
        pass
    if rank > 0:
        returnVal = objective(pop[rank-1],df,str(gen)+"_"+str(rank))
        if returnVal == 0:
            failedBuild += 1
        if returnVal == -1:
            failedRemoval += 1


    if rank > 0:
        buildSend = comm.isend(failedBuild, dest = 0, tag = 4)
        buildSend.wait()
        removalSend = comm.isend(failedRemoval, dest = 0, tag = 5)
        removalSend.wait()

    if rank == 0:
        for i in range(1, size):
            buildRecieve = comm.irecv(source = i, tag = 4)
            failedBuild +=buildRecieve.wait()
            removalRecieve = comm.irecv(source = i, tag = 5)
            failedRemoval += removalRecieve.wait()



    if int(rank) > 0:
        if str(pop[rank-1].getBitList()) in ligandsTested:
            pop[rank-1].setVINA(ligandsTested[str(pop[rank-1])])
            ligandSend = comm.isend(pop[rank-1], dest = 0, tag = 11)
            ligandSend.wait()
        else:
            try:
                score = subprocess.run("./evaluate.sh " + str(gen) + " " + str(rank) + " " + config_path + " " + proteinPDB_path,shell = True, capture_output=True, text=True)
                data = float(score.stdout[0:-1])
                pop[rank-1].setVINA(data)
            except Exception as e:
                pop[rank-1].setVINA(0)

            ligandSend = comm.isend(pop[rank-1], dest = 0, tag = 11)
            ligandSend.wait()
    if rank == 0:
        for i in range(1, size):
            ligandRecieve = comm.irecv(source = i, tag = 11)
            pop[i-1] = ligandRecieve.wait()
    return failedBuild, failedRemoval

# tournament selection
def selection(pop, k=3):
    # first random selection
    random.seed(20)
    selection_ix = randint(len(pop))
    for ix in randint(0, len(pop), k-1):
        # check if better (e.g. perform a tournament)
        if pop[ix].getVINA() < pop[selection_ix].getVINA():
            selection_ix = ix
        
    return pop[selection_ix]

# crossover two parents to create two children
def crossover(p1, p2, r_cross):
    # children are copies of parents by default
    random.seed(20)
    c1 = None
    c2 = None
    newC1, newC2 = p1.getBitList(), p2.getBitList()
    # check for recombination
    if rand() < r_cross:
        # select crossover point that is not on the end of the string
        if len(p1.getBitList()) < len(p2.getBitList()):
            gene = randint(0, len(p1.getBitList()))
        else:
            gene = randint(0, len(p2.getBitList()))
        pt = randint(1, len(p1.getBitList()[gene])-2)
        c1 = p1.getBitList()[gene][:pt] + p2.getBitList()[gene][pt:]
        c2 = p2.getBitList()[gene][:pt] + p1.getBitList()[gene][pt:]
        newC1[gene] = c1
        newC2[gene] = c2
    return [Ligand(newC1), Ligand(newC2)]

# mutation operator
def mutation(ligand, r_mut, r_distant_mut, df):
    #create copy of array (needs to be copy of arrays within array too)
    ligand_copy = []
    for i in ligand.getBitList():
        ligand_copy.append(i.copy())
    
    for i in range(len(ligand_copy)):
        if rand() < r_mut:
            index = random.choices(population = list(df.index.values), weights = list(df['Probability']), cum_weights = None,k=1)
            ligand_copy[i] = [int(x) for x in dec_to_binary(int(index[0]), len(df.index))]
    return ligand_copy

def cleanLigandFragments(ligand, df):
    #This function is to make sure an operation (crossing over or mutation) does not introduce 
    #   a fragment with fewer fragment ends than the fragment it replaces
    #   If this occurs, we remove the fragments that would not be used in BRICS.BUILD
    fragEndsExtra = 0
    last = None
    end = 0
    for i in range(0,len(ligand.getBitList())):
        index = binaryListToInt(ligand.getBitList()[i])
        smilesFrag = df.at[index, "Fragment"]
        fragEnds = smilesFrag.count('*')
        if i == 0:
            fragEndsExtra+=fragEnds
        else:
            fragEndsExtra-=1
            fragEnds-=1
            fragEndsExtra+=fragEnds
        if fragEndsExtra == 0:
            last = i
    
    if last == None:
        return ligand

    return Ligand(ligand.getBitList()[:last+1])
         


# genetic algorithm
def genetic_algorithm(objective, n_bits, n_iter, n_pop, r_cross, r_mut, r_distant_mut, df, config_path, proteinPDB_path):
    random.seed(20)
    fragmentsUsed = set()
    ligandsTested = {}
    pop = []

    medianList = []
    stDevList = []
    bestPerGen = []
    genData = []

    if int(rank) == 0:
        print(len(df))
        for i in range(n_pop):
            pop.append(Ligand([[int(x) for x in dec_to_binary(randint(0,len(df)), len(df))]]))
    # keep track of best solution
    if int(rank) == 0:
        global best
        best = pop[0] 
    totalFailedBuild = []
    failedFragRemoval = []
    # enumerate generations
    for gen in range(n_iter):
        #add all fragments used to set
        if int(rank) == 0:
            for i in pop:
                for j in i.getBitList():
                    binaryString=""
                    for k in j:
                        binaryString+=str(k)
                    fragmentsUsed.add(binaryString)
        
        if rank == 0:
            for i in range(1, size):
                popSend = comm.isend(pop, dest = i, tag = 1)
                popSend.wait()
                ligandsTestedSend = comm.isend(ligandsTested, dest = i, tag = 14)
                ligandsTestedSend.wait()
        if rank > 0:
            popRecieve = comm.irecv(source = 0, tag = 1)
            pop = popRecieve.wait()
            ligandsTested = comm.recv(source = 0, tag = 14)    

        failedBuild, failedRemoval = evaluate(objective, pop, gen, df, ligandsTested, config_path, proteinPDB_path)
        if int(rank) == 0:
            totalFailedBuild.append(failedBuild)
            failedFragRemoval.append(failedRemoval)
            # check for new best solution
            for i in range(len(pop)):   
                for j in pop[i].getBitList():
                    print(binaryListToInt(j), end = " ")
                print(pop[i].getBitList(), pop[i].getVINA(), pop[i].getQED())
                ligandsTested.update({str(pop[i].getBitList()): pop[i].getVINA()})


            genData.append(pop)
            medianList.append(np.mean([ligand.getVINA() for ligand in pop]))
            stDevList.append(np.std([ligand.getVINA() for ligand in pop]))

            for i in range(len(pop)):
                if pop[i].getVINA() < best.getVINA():  
                    best = copy.deepcopy(pop[i])
                    print("Best Orig:",best.getBitList())
            
            # select parents
            selected = [selection(pop) for _ in range(len(pop))]
            
            # create the next generation
            bestParentsList = copy.deepcopy(pop)
            bestScoresList = [ligand.getVINA() for ligand in pop]

            bestScoresList, bestParentsList = (list(t) for t in zip(*sorted(zip(bestScoresList,bestParentsList))))
                
            children = list()
            
            #add mutants of best parents from previous gen
            for i in range(int(len(bestParentsList)*5/8)):
                mutatedParent = Ligand(mutation(bestParentsList[i], r_mut, r_distant_mut, df))
                children.append(mutatedParent)

            if (len(pop) - int(len(bestParentsList)/8)) % 2 == 0:
                end = len(pop) - (len(children) + int(len(bestParentsList)/8))
            else:
                end = len(pop) - (len(children) + int(len(bestParentsList)/8)) + 1

            for i in range(0, end, 2):
                # get selected parents in pairs avoiding identical parents
                if (selected[i].getBitList() == selected[i+1].getBitList()):
                    j = 0
                    while  (selected[i].getBitList() == selected[i+j].getBitList() and i+j < len(selected)-1):
                        j+=1
                    selected[i+j], selected[i+1] = selected[i+1], selected[i+j]
                parent1, parent2 = selected[i], selected[i+1]

                # crossover
                for c in crossover(parent1, parent2, r_cross):
                    # store for next generation
                    children.append(c)

            for child in children:
                if lowerThanMax(child): 
                    index = random.choices(population = list(df.index.values), weights = list(df['Probability']), cum_weights = None,k=1)
                    print(index[0],dec_to_binary(int(index[0]), len(df.index)), dec_to_binary(int(index[0]), len(df.index)), int(dec_to_binary(int(index[0]), len(df.index)),2))
                    newIndexBitList = [int(x) for x in dec_to_binary(int(index[0]), len(df.index))]
                    newBitList = child.getBitList()
                    newBitList.append(newIndexBitList)
                    child = Ligand(newBitList)

            #add best parents from previous gen
            for i in range(n_pop-len(children)):
                children.append(bestParentsList[i])

            print("Len Children", len(children))
            if len(children) > n_pop:
                children.pop()
        
            # replace population
            pop = children
            for i in range(len(pop)):
                pop[i] = cleanLigandFragments(pop[i], df)



            print("Best End Gen " + str(gen)+":", best, "=", best.getVINA())
            bestPerGen.append(best)

    if int(rank) == 0:
        print("Failed Build:", totalFailedBuild)
        print("Failed Removal:", failedFragRemoval)
        print("Unique Ligands Tested:", len(ligandsTested))
        #plotStats(medianList, stDevList, bestPerGen, n_iter, genData)
        writeCSV(genData)

    if int(rank) > 0:
        return [0, 0]
    else:
        return [best, len(fragmentsUsed)]


def syncSendToRanks(pop, distinctTag):
    if rank == 0:
        for i in range(1, size):
            req = comm.isend(pop, dest = i, tag = distinctTag)
            req.wait()
    if rank > 0:
        req = comm.irecv(source = 0, tag = distinctTag)
        pop = req.wait()


def binaryListToInt(list):
    binString = ""
    for i in list:
        binString+=str(i)
    return int(binString, 2)

def writeCSV(genData):
    df = pd.DataFrame({'SMILES':[],
                        'VINA':[],
                        'QED':[],
                        'MULTI':[],
                        'Molecular Weight':[]
                        })
    for gen in genData:
        for ligand in gen:
            if float(ligand.getVINA()) < 0:
                df.loc[len(df.index)] = [ligand.getSMILES(), ligand.getVINA(), ligand.getQED(), ligand.getMULTI(),ligand.getMW()]
    
    df2 = df.sort_values(by = ['VINA','QED'])
    df3 = df2.drop_duplicates(keep='first')

    
    jobID = os.getenv('SLURM_JOB_ID')

    df3.to_pickle("jobResults/results_"+jobID+".pkl")

def plotStats(medianList, stDevList, bestPerGen, n_iter, genData):
    xAxis = []

    vinaAverages = []
    qedAverages = []

    bestVINA = []
    bestQED = []
    bestMULTI = []
    multiAverages = []

    individualVINA = []
    individualQED = []

    for gen in genData:   
        vinaAverages.append(np.mean([ligand.getVINA() for ligand in gen]))
        qedAverages.append(np.mean([ligand.getQED() for ligand in gen]))
        multiAverages.append(np.mean([ligand.getMULTI() for ligand in gen]))
        for ligand in gen:
            individualVINA.append(ligand.getVINA())
            individualQED.append(ligand.getQED())

    for i in bestPerGen:
        bestVINA.append(i.getVINA())
        bestQED.append(i.getQED())
        bestMULTI.append(i.getMULTI())

    for i in range(0,n_iter):
        xAxis.append(i)

    fig, ax = plt.subplots(2,2)

    ax[0,0].scatter(xAxis, vinaAverages)
    ax[0,0].scatter(xAxis, bestVINA)
    ax[0,0].set_title("VINA Scores")

    ax[1,0].scatter(xAxis, qedAverages)
    ax[1,0].scatter(xAxis, bestQED)
    ax[0,0].set_title("QED Scores")

    ax[0,1].scatter(xAxis, multiAverages)
    ax[0,1].scatter(xAxis, bestMULTI)
    ax[0,0].set_title("Multi Scores")

    ax[1,1].scatter(individualVINA, individualQED)
    ax[1,1].set_title("VINA vs QED")
    
    jobID = os.getenv('SLURM_JOB_ID')
    plt.savefig("graphs/gen_stats_"+jobID+".png")

def lowerThanMax(ligand):
    totalSum = 0
    numFragEnds = 0
    for i in ligand.getBitList():
        bitString = ""
        for j in i:
            bitString+=str(j)
        index = int(bitString, 2)
        if index < len(df):
            frag = df.at[index, "Fragment"]
            totalSum += Descriptors.MolWt(Chem.MolFromSmiles(frag))
            for c in frag:
                if c =="*":
                    numFragEnds+=1
    ligand.setMW(totalSum)
    if totalSum < 500 and (numFragEnds-len(ligand.getBitList()) > 0 or len(ligand.getBitList()) == 1):
        return True
    return False

def dec_to_binary(my_int, length):
    number_of_bits = min_bits_to_represent(length)
    """   
    Format a number as binary with leading zeros
    """
    if my_int < length:
        return ("{0:0" + str(number_of_bits) + "b}").format(my_int)
    else:
        return ("{0:0" + str(number_of_bits) + "b}").format(0)

def min_bits_to_represent(size):
    if size == 0:
        return 1
    else:
        return math.floor(math.log2(size))

if __name__ == "__main__":
    global comm
    comm = MPI.COMM_WORLD
    global size
    size = comm.Get_size()
    global rank
    rank = comm.Get_rank()

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--inputcsv", help="Path to fragment CSV from fragmentation pipeline output")
    parser.add_argument("-c","--config", help="Path to Autodock VINA config file specifying search grid for VINA to use")
    parser.add_argument("-p","--protein", help="Path to protein PDBQT file used by Autodock VINA for scoring function")

    args = parser.parse_args()

    df = pd.read_csv(args.inputcsv)
    n_bits = min_bits_to_represent(len(df))
    maxLength = 2**n_bits-1
    df2 = df[:maxLength+1].copy()

    
    #probability function
    for i in range(len(df2)):
        df2.at[i,"MW"] = Descriptors.MolWt(Chem.MolFromSmiles(df2.at[i,"Fragment"]))
        df2.at[i,"Probability"] = ((len(df2)-i)/sum(range(0,len(df2))))*100

    if rank == 0:
        print(df2)

    # define the total iterations
    n_iter = 20
    # define the population size
    n_pop = 40 #MUST BE EVEN NUMBER
    # crossover rate
    r_cross = 0.9
    # mutation rate
    r_mut = 1.0 / float(n_bits) * 5
    r_distant_mut = 1.0 / float(n_bits) #depricated
    # perform the genetic algorithm search
    best, fragmentsUsed = genetic_algorithm(maxVina, n_bits, n_iter, n_pop, r_cross, r_mut, r_distant_mut, df2, args.config, args.protein)
    if int(rank) == 0:
        print('Done!')
        print('f(%s) = %f' % (best, best.getVINA()))
        print("Fragments Used:",fragmentsUsed)
        fragList = []
        for i in best.getBitList():
            binaryString = ''
            for j in i:
                binaryString+=str(j)
            fragList.append(binaryString)
        listBinaryStrings = ' '.join(fragList)


