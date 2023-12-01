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
        frag = targetDict[i]['fragment Key'].split(":")[0]
        if frag not in resultDict:
            details = {}
            details['ligand'] = targetDict[i]['ligand']
            details['domain']=targetDict[i]['domain']
            resultDict[frag] = details
        elif frag in resultDict:
            print("it exists already")
            resultDict[frag]['ligand'].update(targetDict[i]['ligand'])
            resultDict[frag]['domain'].update(targetDict[i]['domain'])
    #print(resultDict)
    df1 = pd.DataFrame.from_dict(resultDict)

    returnDF = df1.transpose()
    returnDF.reset_index(inplace=True)
    returnDF.rename(columns={'index':'fragment'},inplace=True)

    returnDF.to_csv("combinedResidues.csv",index=False)