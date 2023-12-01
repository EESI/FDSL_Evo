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
import parser
import time
import sys
import argparse
import multiprocessing
    
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
            smilesFrag = df.at[index, "SMILES"]
            smilesFragList.append(smilesFrag)
        molFrags = [Chem.MolFromSmiles(f) for f in smilesFragList]
        result = molFrags[0]
        try:
            BRICSObject = BRICS.BRICSBuild(molFrags, maxDepth=1, scrambleReagents=False, onlyCompleteMols=False) #onlyCompleteMols=False to keep frags
            result = list(BRICSObject)[0]
        except:
            pass
        if result != None:
            Chem.SanitizeMol(result)
            return ligand.writeLigand(rank, result, name)
        return 0  

def evaluate(objective, pop, gen, df, ligandsTested, config, proteinPDBQT, proteinPDB):
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


    #prevents ligands over 700 g/mol or ligands which have already been evaluated from running through VINA
    if int(rank) > 0:
        if str(pop[rank-1].getSMILES()) in ligandsTested:
            pop[rank-1].setVINA(ligandsTested[str(pop[rank-1])])
            ligandSend = comm.isend(pop[rank-1], dest = 0, tag = 11)
            ligandSend.wait()
        elif pop[rank-1].getMW() > 700 and pop[rank-1].getPreviousLigand() != None:
            pop[rank-1].setVINA(pop[rank-1].getPreviousLigand().getVINA())
            pop[rank-1].setSMILES(pop[rank-1].getPreviousLigand().getSMILES())
            pop[rank-1].setFinished(True)
            ligandSend = comm.isend(pop[rank-1], dest = 0, tag = 11)
            ligandSend.wait()            
        else:
            try:
                start = time.time() 
                score = subprocess.run("./evaluate.sh " + str(gen) + " " + str(rank) + " " + config + " " + proteinPDBQT + " " + proteinPDB,shell = True, capture_output=True, text=True)
                end = time.time()

                print("Gen:",gen,"| Rank:",rank, "| Time:", end - start)

                pop[rank-1].setGEN(gen)

                data = float(score.stdout[0:-1])
                pop[rank-1].setVINA(data)

            except subprocess.CalledProcessError as exc:
                print("Status : FAIL", exc.returncode, exc.output)
            except Exception as e:
                pop[rank-1].setVINA(0)
                print(e)

            ligandSend = comm.isend(pop[rank-1], dest = 0, tag = 11)
            ligandSend.wait()
    if rank == 0:
        for i in range(1, size):
            ligandRecieve = comm.irecv(source = i, tag = 11)
            pop[i-1] = ligandRecieve.wait()
    return failedBuild, failedRemoval

         

# iterative algorithm
def iterative_algorithm(objective, n_bits, n_iter, n_pop, df, start_index, config, proteinPDBQT, proteinPDB, AADICT, output_directory):
    # initial population of random bitstring
    random.seed(20)
    fragmentsUsed = set()
    ligandsTested = {}
    pop = []

    medianList = []
    stDevList = []
    bestPerGen = []
    genData = []

    TMPDIRNAME = os.getenv('TMPDIR')
    subprocess.run("mkdir " + TMPDIRNAME + "/pdbForPLIP > /dev/null", shell = True)
    subprocess.run("mkdir " + TMPDIRNAME + "/pdbForPLIP/"+str(rank) + " > /dev/null", shell = True)
    subprocess.run("cp plip_v2.2.0.simg " +  TMPDIRNAME + "/pdbForPLIP/"+str(rank) + " > /dev/null", shell = True)


    if int(rank) == 0:
        print(len(df))
        for i in range(n_pop):
            pop.append(Ligand([[int(x) for x in dec_to_binary(start_index,len(df))]], df.at[start_index,"SMILES"]))
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
            pop = comm.recv(source = 0, tag = 1)
            ligandsTested = comm.recv(source = 0, tag = 14)
    

        failedBuild, failedRemoval = evaluate(objective, pop, gen, df, ligandsTested, config, proteinPDBQT, proteinPDB)

        #determines whether ligands are valid
        if rank > 0:
            print("STARTING VALIDATION FOR", rank, pop[rank-1].getSMILES())
            validation = validate_smiles_for_output(pop[rank-1].getSMILES())
            send = comm.isend(validation, dest = 0, tag = 6) 
            send.wait()

        else:
            print("Rank 0 STARTING VALIDATION")
            validLigandList = []
            for i in range(1,size):
                print("RECIEVING VALIDATION FOR", i)
                recieve = comm.irecv(source = i, tag = 6)
                results = recieve.wait()
                if results and not pop[i-1].isFinished():
                    pop[i-1].setTASK(os.getenv("SLURM_ARRAY_TASK_ID"))
                    pop[i-1].setRANK(i)
                    pop[i-1].setGEN(gen)
                    validLigandList.append(pop[i-1])
                else:
                    print("FAILED VALIDATION")
                
            genData.append(validLigandList)
            print("FINISHED VALIDATION")

            

        if int(rank) == 0:
            totalFailedBuild.append(failedBuild)
            failedFragRemoval.append(failedRemoval)
            # check for new best solution
            
            for i in range(len(pop)):   
                
                pop[i].setMW(df)
            
                for j in pop[i].getBitList():
                    print(binaryListToInt(j), end = " ")
                print(pop[i].getSMILES(), pop[i].getVINA(), pop[i].getQED(), pop[i].getMW())

                ligandsTested.update({pop[i].getSMILES(): pop[i].getVINA()})

            medianList.append(np.mean([ligand.getMULTI() for ligand in pop]))
            stDevList.append(np.std([ligand.getMULTI() for ligand in pop]))

            for i in range(len(pop)):
                if pop[i].getMULTI() < best.getMULTI():  
                    best = copy.deepcopy(pop[i])
                    print("Best Orig:",best.getSMILES())
            
            #GENERATE NEXT GENERATION
            errors = ""
            
            #prints all ligands in a generation and VINA scores
            for i in range(len(pop)):
                try:
                    if pop[i].getPreviousLigand() == None or pop[i].getVINA() < pop[i].getPreviousLigand().getVINA():
                        if pop[i].getPreviousLigand() != None:
                            print("CURRENT BETTER | Gen:", gen, "| Rank:", i+1, "| Current:", pop[i].getVINA(),"| Previous:", pop[i].getPreviousLigand().getVINA())
                        else:
                            print("NONE IN PREVIOUS | Gen:", gen, "| Rank:", i+1)
                        pop[i].setPreviousLigand(copy.deepcopy(pop[i]))
                    elif pop[i].getPreviousLigand() != None:
                        if pop[i].getVINA() == pop[i].getPreviousLigand().getVINA():
                            print("ADDITION FAILED | Gen:", gen, "| Rank:", i+1, "| Current:", pop[i].getVINA(),"| Previous:", pop[i].getPreviousLigand().getVINA())
                        else:
                            print("CURRENT WORSE | Gen:", gen, "| Rank:", i+1, "| Current:", pop[i].getVINA(),"| Previous:", pop[i].getPreviousLigand().getVINA())
                    else:
                        print("NOT REACHED | Gen:", gen, "| Rank:", i+1, "| Current:", pop[i].getVINA())
                except Exception as e:
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    print(exc_type, fname, exc_tb.tb_lineno)

            print("Best End Gen " + str(gen)+":", best, "=", best.getVINA())
            bestPerGen.append(best)

        try:
            pop = comm.bcast(pop, root = 0)
        except Exception as e:
            print("ERROR AT RANK:", rank, e, "POP:",pop)

        if rank > 0:
            if pop[rank-1].getMW() < 700:
                try:
                    output = subprocess.run("./addFragment.sh " + str(gen) + " " + str(rank) + " " + str(pop[rank-1].getPreviousLigand().getGEN()) + " " + proteinPDB + " " + AADICT,shell = True, capture_output=True, text=True)
                    print(output.stdout)
                    pop[rank-1].reset(currentMol = Chem.MolFromPDBFile(os.getenv('TMPDIR') + "/pdbForPLIP/" + str(rank) + "/fragment_addition_result_" + str(gen) + "_" + str(rank) + ".pdb"), df = df)          
                    print("Rank", rank, "NEW LIGAND:", pop[rank-1].getSMILES())
                    print("Rank", rank, "PREVIOUS LIGAND:", pop[rank-1].getPreviousLigand().getSMILES())
                except Exception as e:
                    print("Rank",rank,"ERROR:",e)
            else:
                print("Rank", rank, "TOO LARGE FOR ADDITION:", pop[rank-1].getSMILES())
            
            send = comm.isend(pop[rank-1], dest = 0, tag = 82) 
            send.wait()
        else:
            validLigands = []
            for i in range(1,size):
                recieve = comm.irecv(source = i, tag = 82)
                pop[i-1] = recieve.wait()

    if int(rank) == 0:
        print("Failed Build:", totalFailedBuild)
        print("Failed Removal:", failedFragRemoval)
        print("Unique Ligands Tested:", len(ligandsTested))
        #plotStats(medianList, stDevList, bestPerGen, n_iter, genData)
        writeCSV(genData, output_directory)

    if int(rank) > 0:
        return [0, 0]
    else:
        return [best, len(fragmentsUsed)]

def run_function_with_timeout(func, timeout, *args, **kwargs):
    # Create a new process
    process = multiprocessing.Process(target=func, args=args, kwargs=kwargs)

    # Start the process
    process.start()

    # Wait for 'timeout' seconds or until the process terminates
    process.join(timeout)

    # If the process is still alive, terminate it and return False
    if process.is_alive():
        process.terminate()
        process.join()
        return False

    # If the process completed within the time limit, return its result
    return process.exitcode == 0

def validate_smiles_for_output(smiles):
    try:
        mol2D = Chem.MolFromSmiles(smiles)
        molFromSmilesWithH = Chem.AddHs(mol2D)
        AllChem.EmbedMolecule(molFromSmilesWithH, randomSeed=739)   #add standardized seed
        AllChem.MMFFOptimizeMolecule(molFromSmilesWithH)
        embeddedMolecule = Chem.RemoveHs(molFromSmilesWithH)
        return True
    except Exception as e:
        return False

def writeCSV(genData, output_directory):
    df = pd.DataFrame({'SMILES':[],
                        'VINA':[],
                        'QED':[],
                        'MULTI':[],
                        'Molecular Weight':[],
                        'Task':[],
                        'Generation':[],
                        'Rank':[]
                        })
    for gen in genData:
        for ligand in gen:
            if float(ligand.getVINA()) < 0:
                df.loc[len(df.index)] = [ligand.getSMILES(), ligand.getVINA(), ligand.getQED(), ligand.getMULTI(),ligand.getMW(),ligand.getTASK(), ligand.getGEN(), ligand.getRANK()]

    df.sort_values(by = ['VINA','QED'])

    df2 = df.drop_duplicates(keep='first')

    jobID = os.getenv('SLURM_JOB_ID')

    df2.to_pickle(output_directory + "results_"+jobID+".pkl")


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

def sameAAInteractions(ligand, newIndex, df):
    if ligand.getAAInteractions() != None:
        for AA in ligand.getAAInteractions():
            if AA in df.at[newIndex, "domain"]:
                return True
    return False

def validAddition(ligand, df):
    ligand.setMW(df)
    if ligand.getMW() < 500 and (ligand.getFragmentEnds()-len(ligand.getBitList()) > 0 or len(ligand.getBitList()) == 1):
        return True
    return False

def dec_to_binary(my_int, length):
    """   
    Format a number as binary with leading zeros
    """
    if my_int < length:
        return "{0:018b}".format(my_int)
    else:
        return "000000000000000000"



if __name__ == "__main__":
    global comm
    comm = MPI.COMM_WORLD
    global size
    size = comm.Get_size()
    global rank
    rank = comm.Get_rank()

    argparser = argparse.ArgumentParser()

    argparser.add_argument("-i", "--index", help="Index in ligand list to use as starting ligand")
    argparser.add_argument("-s", "--starterligands", help = "CSV used as starting ligands")
    argparser.add_argument("-c","--config", help="Path to Autodock VINA config file specifying search grid for VINA to use")
    argparser.add_argument("-p","--proteinPDBQT", help="Path to protein PDBQT file used by Autodock VINA for scoring function")
    argparser.add_argument("-x","--proteinPDB", help="Path to protein PDB file used for PLIP")
    argparser.add_argument("-a","--aadictionary", help="Path to AADictionary.pickle outputted from generateFragmentBank.sh")
    argparser.add_argument("-o","--outputdirectory", help="Path to output directory where pickle files should be stored")

    args = argparser.parse_args()
    df = pd.read_csv(args.starterligands)

    df2 = df[df['Molecular Weight'] < 550] 

    df3 = df2.reset_index()

    if rank == 0:
        print(df3)

    # define the total iterations
    n_iter = 10
    # bits
    n_bits = 18
    # define the population size
    n_pop = 40
    start_index = int(args.index)

    best, fragmentsUsed = iterative_algorithm(maxVina, n_bits, n_iter, n_pop, df3, start_index, args.config, args.proteinPDBQT, args.proteinPDB, args.aadictionary, args.outputdirectory)
    if int(rank) == 0:
        print('Done!')
        print('f(%s) = %f' % (best, best.getVINA()))
        print("Fragments Used:",fragmentsUsed)
