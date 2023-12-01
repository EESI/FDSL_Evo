#!usr/bin/python

import os
import sys



if __name__=="__main__":

    if (len(sys.argv)!=5):
        print('Need exactly 4 input: the dir where ligands is located (first arg), and the main protien file, as well as the tmp directory, and the ID of the job array')
    else:
        ligand_dir= sys.argv[1]
        protein = sys.argv[2]
        TMP=sys.argv[3]
        ID=sys.argv[4]
        for filename in os.listdir(ligand_dir):
            print('this file is being processed', file = sys.stdout)
            print(filename, file = sys.stdout)
            target = os.path.join(ligand_dir,filename)
            command = "python "+ TMP+"/merge-it.py -ir " + protein + " -if " + target
            os.system(command)
        os.system("mv "+TMP+"/"+ID+"RelA-dehyd-* "+ TMP+"/Combined/"+ID)
        print("mv "+TMP+"/"+ID+"RelA-dehyd-* "+ TMP+"/Combined/"+ID) 
