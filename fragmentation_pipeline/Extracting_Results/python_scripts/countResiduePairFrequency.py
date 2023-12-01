

import os
import sys
import pandas as pd
import pickle

if __name__=="__main__":
    if len(sys.argv) < 2:
        sys.exit()
    filename=sys.argv[1]
    # filename = whatever the main.pickle is
    targetDict=pickle.load( open(filename,'rb'))
    theDict={}
    for ligand in targetDict:
        FragDict = targetDict[ligand]['FragToResDict']
        originalFragDict = targetDict[ligand]['originalFragments']
        for frag in FragDict:
            fragKey = frag + ": " + ",".join(sorted(FragDict[frag]))
            originalFrag=originalFragDict[frag]
            originalFragKey = originalFrag + ": " + ",".join(sorted(FragDict[frag]))
            if originalFragKey not in theDict:
                details = {}
                details['ligand'] = set()
                details['ligand'].add(ligand)
                details['fragment']=frag
                details['domain']=FragDict[frag]
                theDict[originalFragKey]=details
            elif originalFragKey in theDict:
                print("it exists already")
                theDict[originalFragKey]['ligand'].add(ligand)
    theDF = pd.DataFrame.from_dict(theDict)

    newDF=theDF.transpose()
    #theDF.drop('placeHolder',inplace=True,axis=1)
    newDF.reset_index(inplace=True)
    newDF.rename(columns={'index':'fragment Key'},inplace=True)

    print(newDF)
    
    newDF.to_pickle("residueFrequency.pickle")
    newDF.to_csv("residueFrequency.csv",index=False)
                
