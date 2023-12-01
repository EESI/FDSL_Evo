#!/usr/bin/python


import os
import sys


def processExcessoryFile(needToDestroy,ligandDir,filename):
    '''
    removed all the unnecessary files to prevent them from confusing the program
    requires 3 inputs(BAD PROGRAMMING): needToDestroy's file name + the ligand's directory + the filename itself

    '''
    os.system('rm ' + needToDestroy)
    os.system('rm ' + ligandDir + 'plipfixed* ')
    os.system('mv ' + ligandDir + 'report.xml ' + ligandDir + "xmls/" + filename + '.xml')


def runSingularity(ligandDir):
    '''
    Runs the singularity module that would then run the command in os.
    input-specs: requires a ligand directory
    output- all output resides in the xml directory. 
    '''
    commandString = 'singularity run --bind ' + ligandDir + ':/mnt' + ' plip_v2.2.0.simg'
    os.system('mkdir ' + ligandDir + 'xmls')
    for filename in os.listdir(ligandDir):
        newCommand = commandString + " " + "-f /mnt/" + filename + " -o /mnt -x "
        print(newCommand)
        os.system(newCommand)
        needToDestroy = ligandDir
        splittedName = filename.split("-")
        needToDestroy += "-".join(splittedName[0:3])
        needToDestroy += '-prep_protonated.pdb'
        processExcessoryFile(needToDestroy, ligandDir,filename)


if __name__=="__main__":
    if (len(sys.argv)) != 2:
        print("need only 1 argument")
    else:
        runSingularity(sys.argv[1])
    # os.system('singularity run --bind .:/mnt plip_v2.2.0.simg -f /mnt/1.pdb -o /mnt -x ')


