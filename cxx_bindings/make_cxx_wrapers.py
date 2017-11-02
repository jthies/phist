#!/usr/bin/python
import re       # Regular expression library
import os
import time
import csv      # To write .csv files
import sys      # To use command line arguments
import platform  # To obtain platform informations
from collections import OrderedDict


library_path="../src/"
for directory in os.walk(library_path):
    for c_file_name in directory[2]:
        if (c_file_name.find("_decl") == -1):
            continue

        whole_file_original = open(directory[0] + "/" + c_file_name, 'rU').read()
        whole_file = whole_file_original

        subroutines=re.findall("SUBR\(.*\)\(.*\)\;", whole_file)
        print subroutines
