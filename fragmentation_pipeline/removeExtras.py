import pandas as pd
import os
import sys




LetterToDict = {'A': 2,'B': 3,'C': 5, 'D': 7, 'Unidentified Domain': 11}


def convertSetToInt(theSet):
    '''
        This converts the subdomain's information (which is a set of Domain represented by letters) to an integer value
        How it works: It maps each letter to a prime number, and then multiply the numbers if the letter is in the set of domain
    '''
    num=1
    for i in theSet:
         if i in LetterToDict:
            num *= LetterToDict[i]   
    return num

if __name__=="__main__":
    theDF = pd.read_csv("3bikTestNoExtra.csv")
    #theDF["domain_int"] = 1
    #theDF["domain"] = theDF["domain"].apply(lambda x: eval(x))
    #theDF["domain_int"] = theDF["domain"].apply(lambda x: convertSetToInt(x))

    #theDF.drop(columns=["Fragment","domain"],inplace=True)
    theDF.drop(columns=["Fragment"],inplace=True)


    print(theDF)
    theDF.to_csv("resultWithDomainInfo.csv",index=False)
