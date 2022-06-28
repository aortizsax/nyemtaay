#! /usr/bin/env python
# -*- coding: utf-8 -*-
# two functions for parsing a list of fasta files and a meta data

###########################################################################
## Reader Implementation

import dendropy
import numpy as np
import pandas as pd


def read_population_data(filename):
    print(pd.read_csv(filename[0], header=0, index_col=0))
    return

def polymorphism_read_csv(filename,header):
    print(dendropy.DnaCharacterMatrix.get(
                    path=filename,
                    schema="fasta",
                    label=None,
                    taxon_namespace=None,
                    matrix_offset=None,
                    ignore_unrecognized_keyword_arguments=False,
                    ))   
    return

def read_fasta_files(filenames):
    
    fasta_dict = {}
    
    for file in filenames:
        with open(file) as f:
            lines = f.readlines()
            taxa = lines[0].strip()
            fasta_dict[taxa[1:]] = lines[1]
    length_sequnces = len(lines[1])
    
    dna = dendropy.DnaCharacterMatrix.from_dict(fasta_dict)
    
    print(dna.taxon_state_sets_map(char_indices=range(length_sequnces), 
          gaps_as_missing=False, 
          gap_state=None, 
          no_data_state=None))

    # with open(filenames[0]) as f:
        # lines = f.readlines()
        # nt = lines[1].strip()
        # width_matrix = len(nt) #

    # # make arrary  null
    # data = np.zeros([height_matrix, width_matrix], dtype="int64")
    # options = {"A": 0, "G": 1, "T": 2, "C": 3}

    # for i, file in enumerate(filenames):
        # print('new file')
        
        # with open(file) as f:
            # lines = f.readlines()
            # cc = 0
            # for j,nt in enumerate(lines[1]):
                # char = options[nt]
                # data[i, j] = char
                # cc += 1
    # print(data)
    
    return