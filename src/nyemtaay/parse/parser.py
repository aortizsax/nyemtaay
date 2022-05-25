#! /usr/bin/env python
# -*- coding: utf-8 -*-
# two functions for parsing a list of fasta files and a meta data

###########################################################################
## Reader Implementation

import dendropy
import numpy as np
import pandas as pd


def read_metadata(filename=""):
    print(pd.read_csv(filename, header=0, index_col=0))
    return


def read_fasta_files(filenames=[]):

    li = len(filenames)
    lj = 1

    # make arrary  null
    data = np.zeros([li, lj], dtype="int64")
    options = {"A": 0, "G": 1, "T": 2, "C": 3}

    for i, file in ennumerate(filenames):
        with open(file) as f:
            lines = f.readlines()
            cc = 0
            for j,nt in enumerate(lines[1]):
                char = options[nt]
                data[i, j] = char
                cc += 1

    print(data)
    return


pass
