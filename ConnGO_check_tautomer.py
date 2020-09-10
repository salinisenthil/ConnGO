import os, sys
import linecache
import numpy as np


file_name = sys.argv[1]
dir_path = sys.argv[3]
or_smi = sys.argv[2]
new_fname = dir_path + "/modified_tautomers.smi"
or_path = dir_path + "/" + or_smi

for line in open(or_path,"r"):
    if file_name in line:
        with open("temp_file","w") as tf:
            tf.write(line.strip())

for line in open("temp_file", "r"):
        if "[NH3+]" in line:
            rp1 = line.replace("[NH3+]","N")
            rp1 = rp1.replace("[O-]", "O")
            with open(new_fname, "a") as nf:
                nf.write(rp1 + "\n") 
        if "[NH2+]" in line:
            rp1 = line.replace("[NH2+]","N")
            rp1 = rp1.replace("[O-]", "O")
            with open(new_fname, "a") as nf:
                nf.write(rp1 + "\n") 
        if "[NH+]" in line:
            rp1 = line.replace("[NH+]","N")
            rp1 = rp1.replace("[O-]", "O")
            with open(new_fname, "a") as nf:
                nf.write(rp1 + "\n") 

