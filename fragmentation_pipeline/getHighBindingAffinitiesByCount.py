import pandas as pd
import sys
import os

if __name__ == "__main__":
    if len(sys.argv)!=2:
        print("need csv file name to alter")
        sys.exit()
    else:
        csvData = pd.read_csv(sys.argv[1])
        print("unsorted")
        #print(csvData)
        csvData['Mean'] = csvData['Mean'].astype(float)
        csvSorted = csvData.sort_values(["Mean"],axis=0,ascending=[True],inplace=False)
        print("sorted")
        #print(csvSorted)
        top25Percent = csvSorted.head(int(len(csvSorted)*(25/100)))
        print("top 25%")
        #print(top25Percent)
        top10Percent = csvSorted.head(int(len(csvSorted)*(10/100)))
        print("top 10%")
        #print(top10Percent)
        top25PercentSorted = top25Percent.sort_values(["Count"],axis=0,ascending=[False],inplace=False)
        print("top 25 by count")
        #print(top25PercentSorted)        
        top10PercentSorted = top10Percent.sort_values(["Count"],axis=0,ascending=[False],inplace=False)
        print("top 10 by count")
        #print(top10PercentSorted) 
        top25PercentSorted.to_csv("3F4MsortedResults/" + sys.argv[1].split(".")[0] + "Top25Percent.csv",index=False)
        top10PercentSorted.to_csv("3F4MsortedResults/" + sys.argv[1].split(".")[0] + "Top10Percent.csv",index=False)
