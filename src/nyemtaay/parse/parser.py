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
## DISCLAIMED. IN NO EVENT SHALL ADRIAN ORTIZ-VELEZ BE LIABLE FOR ANY DIRECT,
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


def read_metadata(filename, hd):
    population_structure = pd.read_csv(filename, header=int(hd), index_col=0)
    # print(population_structure)
    return population_structure


def read_fasta_files(filenames):
    if len(filenames) > 1:
        number_of_files = len(filenames)

        with open(filenames[0]) as f:
            lines = f.readlines()
            nt = lines[1].strip()
            length_sequence = len(nt)

        # make arrary  null
        data = np.zeros([number_of_files, length_sequence], dtype="int64")
        alphabet = {"A": 0, "G": 1, "T": 2, "C": 3, "R": 4}

        for i, file in enumerate(filenames):
            with open(file) as f:
                lines = f.readlines()
                for j, nt in enumerate(lines[1]):
                    char = alphabet[nt]
                    data[i, j] = char

        return data

    else:
        with open(filenames[0]) as f:
            lines = f.readlines()
            nt = lines[1].strip()
            labels = []

            for line in lines:
                if line.startswith(">"):
                    labels.append(line)

            label_count = -1
            alignment = [""] * len(labels)
            for line in lines:
                if line.startswith(">"):
                    label_count += 1
                else:
                    alignment[label_count] += line.strip()

            length_sequence = len(alignment[0])
            number_sequences = len(labels)

            # make arrary  null
            data = np.zeros([number_sequences, length_sequence], dtype="int64")
            alphabet = {
                "A": 0,
                "G": 1,
                "T": 2,
                "C": 3,
                "R": 4,
                "Y": 4,
                "-": 4,
                "K": 4,
                "H": 4,
                "W": 4,
                "N": 4,
                "S": 4,
                "M": 4,
                "B": 4,
                "V": 4,
                "D": 4,
            }

            for i, label in enumerate(labels):
                label = label.strip()[1:]

                alignment[i] = alignment[i].strip()
                for j, SNP in enumerate(alignment[i]):
                    data[i, j] = alphabet[SNP]

        #        print('length',len(alignment[0]),'num sequenes',len(labels))
        return data
        
def convert_to_fasta():
    #samtools fastq -0 /dev/null in_name.bam > all_reads.fq
    #pysam.fasta('-0', "/dev/null", "in_name.bam")
    #how to > 
    
    #f = open("in_name.fasta", "w+")
    #f.write(pysam.fasta('-0', "/dev/null", "in_name.bam"))
    #f.close()

    ##open and read the file after the appending:
    #f = open("in_name.fasta", "r")
    #print(f.read())
    return

def to_dataframe(sequences, data):
    index = data.index
    length = len(sequences[0])
    sequence_dataframe = pd.DataFrame(
        sequences, columns=range(1, length + 1), index=index
    )
    return sequence_dataframe


pass
