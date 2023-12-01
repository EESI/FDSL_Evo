import pandas as pd 
import xmltodict 
import ast 
import os
from json import loads, dumps
from collections import OrderedDict
import re
from ast import literal_eval
from  itertools import chain 
import sys

if __name__ == "__main__":

    xml_path_dir = sys.argv[1]
    xml_save_dir = sys.argv[2]
    for filename in os.listdir(xml_path_dir):
        if filename.endswith(".xml"):
            xml_path=os.path.join(xml_path_dir,filename)
            xml_string = open(xml_path).read()
            binding = re.search("<bindingsite id=\"2\"", xml_string)
            if binding == None:
                print(filename)
            else:
                binding_start=binding.start()
                bindingFile = xml_string[binding_start:-1]
                ligandName=filename.split("-")[2]
                newFileName = ligandName+".xml"
                newFileNamePath=os.path.join(xml_save_dir,newFileName)

                with open(newFileNamePath,"w") as f:
                    for line in bindingFile:
                        f.write(line)
        

