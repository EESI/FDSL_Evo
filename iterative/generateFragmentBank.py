import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import pickle
from rdkit.Chem import PandasTools
from rdkit.Chem import AllChem
from rdkit.Chem import SDWriter
import time
import argparse

def filterFragments(csvPath, maxWeight, maxPoolLength):
    df = pd.read_csv(csvPath)

    df = df.sort_values(by='Median', ascending=True)

    fullNameAAs = []
    print(df)
    # Create a new dictionary with keys as numbers and associated lists of letters
    new_data = {}

    for _, row in df.iterrows():
        if Descriptors.MolWt(Chem.MolFromSmiles(row['fragment'])) < maxWeight:
            AAs = row['domain']
            fragment = row['fragment']
            for AA in AAs.replace(" ", "").split(","):
                if AA[3:] not in new_data:
                    new_data.update({AA[3:]: set()})
                    new_data[AA[3:]].add(fragment)
                    fullNameAAs.append(AA)
                elif len(new_data[AA[3:]]) < maxPoolLength:
                    new_data[AA[3:]].add(fragment)

    for key in new_data:
        print(key, len(new_data[key]))
        
    print(new_data)

    return new_data, fullNameAAs


if __name__ == "__main__":
    pd.set_option('display.max_columns', None)

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputcsv", help="Path to NoExtra fragment CSV from fragmentation pipeline output")
    args = parser.parse_args()

    inputCSVPath = args.inputcsv
    maxFragmentSize = 200
    maxPoolLength = 1000
    AADict, fullNameAAs= filterFragments(inputCSVPath, maxFragmentSize, maxPoolLength)
    with open('AADictionary.pickle', 'wb') as handle:
        pickle.dump(AADict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    