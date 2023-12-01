import os
import pandas as pd
import sys

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit()
    filename=sys.argv[1]
    
    targetDict = pd.read_pickle(filename)
    targetDict = targetDict.transpose()
    #print(targetDict)

    #print(targetDict)
    resultDict = {}
    for i in targetDict:
        fragOne = targetDict[i]['frag 1']
        if fragOne not in resultDict:
            details = {}
            details['ligand'] = targetDict[i]['ligand']
            details['frag 2'] = set()
            details['frag 2'].add(targetDict[i]['frag 2'])
            resultDict[fragOne] = details
        elif fragOne in resultDict:
            print("it exists already")
            resultDict[fragOne]['ligand'].update(targetDict[i]['ligand'])
            resultDict[fragOne]['frag 2'].add(targetDict[i]['frag 2'])
    #print(resultDict)
    df1 = pd.DataFrame.from_dict(resultDict)

    returnDF = df1.transpose()
    returnDF.reset_index(inplace=True)
    returnDF.rename(columns={'index':'frag 1'},inplace=True)

    returnDF.to_csv("combinedOrigFragOne.csv",index=False)