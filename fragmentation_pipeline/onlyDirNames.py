

import os
import sys


if __name__ == "__main__":
    """
        Input: first argument is the file "listOfNames"
    """
    theFile = sys.argv[1]
    f2 = open("dirNames.txt","a")

    with open(theFile) as f:
        for line in f:
            arr=line.split("/")
            f2.write(arr[2])
    f2.close()

