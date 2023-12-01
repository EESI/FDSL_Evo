import os
import pickle
import sys
import pandas as pd
import numpy as np
from scipy.stats import ttest_1samp
from scipy import stats
from statsmodels.stats import weightstats as stests
import scipy.stats
from scipy.stats import wasserstein_distance
from scipy.stats import ttest_ind
import itertools

if len(sys.argv) < 5:
    sys.exit()
combineCSV=sys.argv[1]
origFragPairFrequencyCSV=sys.argv[2]
parentDirectory=sys.argv[3]
receptorName=sys.argv[4]

df=pd.read_csv(combineCSV)
fragmentResiduePairDF=pd.read_csv(origFragPairFrequencyCSV)


#Setting up
fragmentResiduePairDF['ligandSet']=fragmentResiduePairDF['ligand'].apply(lambda x: eval(x))
df['FragmentSet']=df['Fragments'].apply(lambda x: eval(x))
data=[df['Ligand'],df['VINA']]
#len(set(itertools.chain.from_iterable(df['FragmentSet'])))
headers=['Ligand','Binding Affinity']
df3 = pd.concat(data, axis=1, keys=headers)
df10=df3.dropna()

#For Fragment Residue Pair
#meanOfEntire = df3['Binding Affinity'].sum()/585281 
#print(meanOfEntire)
targets= fragmentResiduePairDF['fragment Key'].to_list()
theDict={}
entireMedian=df10.median()[0]

for target in targets:
    info={}
    theIndex=fragmentResiduePairDF.loc[fragmentResiduePairDF['fragment Key'] == target].index
    df5=df10[df10['Ligand'].isin(fragmentResiduePairDF.loc[theIndex[0],'ligandSet'])]
    #ztest ,pval = stests.ztest(df5['Binding Affinity'], df3['Binding Affinity'], value=0)
    info['frag 1']= target.split(":")[0]
    info['frag 2']=target.split(":")[1] #added
    #info['domain']=fragmentResiduePairDF.loc[theIndex[0],'domain']
    try:

        mean = df5.mean()[0]
        print("this is the mean",mean)
        mode = df5['Binding Affinity'].mode()[0]
    #mode=5
        median = df5.median()[0]
        info['Mean'] = mean
        info['Mode'] = mode
        info['Median']= median
        info['Diff_Median'] = median - entireMedian
        info['Count'] = len(df5.index)    
        theDict[target] = info
    except:
        print("This is DF5")
        print(df5)
dataDF=pd.DataFrame.from_dict(theDict,orient='index')
dataDF.reset_index(inplace=True)
dataDF.rename(columns={'index':'Key'},inplace=True) #changed Fragment to Key
#dataDF.drop('Fragment',axis=1,inplace=True)

os.chdir(parentDirectory)
dataDF.to_csv(receptorName + 'GenericToSpecificResults.csv',index=False)
