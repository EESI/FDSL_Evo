


import os 
import sys
import pickle
import glob
import pandas as pd


if __name__ == "__main__":
    mainPickleFolder = sys.argv[1]
    os.chdir(mainPickleFolder)
    # the argument should be the absolute path towards the "main_pickle" directory
    extension = 'pickle'
    all_filenames = [i for i in glob.glob('*.{}'.format(extension))]
    mainDict={}
    for pick in all_filenames:
        try:
            targetDict=pickle.load( open(pick,'rb'))
            for ligand in targetDict:
                mainDict[ligand]=targetDict[ligand]
        except:
            print("file empty")
    pickle.dump(mainDict,open('combine_pickle.pickle','wb')) 
    df=pd.DataFrame.from_dict(mainDict)
    df2=df.transpose()
    df2.reset_index(inplace=True)
    df2.rename(columns={'index':'Ligand'},inplace=True)
    print(df2)
    df2.sort_values(by='VINA',inplace=True,ascending=False)
    df2.to_csv('combine.csv',index=False)

