#! /usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
## Copyright (c) 2022 Adrian Ortiz-Velez.
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.
##     * Redistributions in binary form must reproduce the above copyright
##       notice, this list of conditions and the following disclaimer in the
##       documentation and/or other materials provided with the distribution.
##     * The names of its contributors may not be used to endorse or promote
##       products derived from this software without specific prior written
##       permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
## ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
## WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL JEET SUKUMARAN BE LIABLE FOR ANY DIRECT,
## INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
## BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
## LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
## OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
## ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
##############################################################################
#
# two functions for parsing a list of fasta files and a meta data
#
###########################################################################
## Reader Implementation

import dendropy
import numpy as np
import pandas as pd


def read_metadata(filename):
    population_structure = pd.read_csv(filename, header=0, index_col=0)
    print(population_structure)
    return population_structure


def read_fasta_files(filenames):
    
    print(len(filenames))
    
    if len(filenames) > 1:
        number_of_files = len(filenames)
        
        with open(filenames[0]) as f:
            lines = f.readlines()
            nt = lines[1].strip()
            length_sequence = len(nt)        
            
        # make arrary  null
        data = np.zeros([number_of_files, length_sequence], dtype="int64")
        alphabet = {"A": 0, "G": 1, "T": 2, "C": 3}

        for i, file in enumerate(filenames):
            with open(file) as f:
                lines = f.readlines()
                for j,nt in enumerate(lines[1]):
                    char = alphabet[nt]
                    data[i, j] = char

        print(data)
        return data
    
    else:
        with open(filenames[0]) as f:
            lines = f.readlines()
            nt = lines[1].strip()
            length_sequence = len(nt)
            number_sequences = int(len(lines)/2)
            
            # make arrary  null
            data = np.zeros([number_sequences, length_sequence], dtype="int64")
            alphabet = {"A": 0, "G": 1, "T": 2, "C": 3}
            
            labels = lines[::2]
            alignment = lines[1::2]
            print(labels)
            print(alignment)
            
            for i, label in enumerate(labels):
                label = label.strip()[1:]
                print(i,label)
                alignment[i] = alignment[i].strip()
                for j, SNP in enumerate(alignment[i]):
                    data[i,j] = alphabet[SNP]
                    
        print('length',length_sequence,'num sequenes',number_sequences)
        return data
pass














