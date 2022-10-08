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

# Libraries 
import numpy as np 
import matplotlib.pyplot as plt
import dendropy
import pandas as pd
from nyemtaay.mathlib import mathfn

def segregating_sites(sequence_dataframe,data):
    print("identifying segratating sites")
    
    seg_sites = []
    for position in sequence_dataframe:
        frequncies = sequence_dataframe[position].value_counts(normalize=True)
        if len(frequncies) > 1:
            seg_sites.append(position)
            
    print(seg_sites)
    return seg_sites

def number_segregating_sites(sequence_dataframe,data):
    print("calculating number segratating sites")
    seg_sites = segregating_sites(sequence_dataframe,data)
    print("number segrating sites:",len(seg_sites))
    return len(seg_sites)
    
def number_pairwise_differences(sequence_dataframe):

    count_pairwise_differences = 0
    
    
    for i, (index, sequence) in enumerate(sequence_dataframe[:-1].iterrows()): 
        #print(i,sequence)
        difference = abs(sequence - sequence_dataframe[i+1:])
        #print(difference)
        normalized_diff = (difference / difference)
        
        count_pairwise_differences += normalized_diff.sum().sum()
    return count_pairwise_differences
    
def nucleotide_diversity(sequence_dataframe,data):
    print("calculating nucleotide diversity(pi)")
    

    number_sequences = len(sequence_dataframe)
    number_comparisons = mathfn.nchoose2(number_sequences)
    length = len(sequence_dataframe.iloc[0])
    
    nucleotide_div = number_pairwise_differences(sequence_dataframe)
    nucleotide_div /= number_comparisons
    nucleotide_div /= length
    
    print("nucdiv",nucleotide_div)
    return nucleotide_div
    
    




def frequency_spectrum(sequence_dataframe,data):
    print("calculating sfs (add per sequnce length)")
    length = len(sequence_dataframe.iloc[0]) 
    
    site_frequency_spectrum = [0] * int(length/2)
    for position in sequence_dataframe:
        site_counts = sequence_dataframe[position].value_counts().tolist()
        if len(site_counts) > 1:
            site_frequency_spectrum[site_counts[1]] += 1
            
    print(site_frequency_spectrum)
    
    y_pos = np.arange(len(site_frequency_spectrum))

    plot_length = int(0.05*len(site_frequency_spectrum))
    
    plt.bar(y_pos[:plot_length],site_frequency_spectrum[:plot_length])
    plt.show()
    
    return site_frequency_spectrum

def frequency(sequence_dataframe,data):
    print("calculating frquencies")
    print(sequence_dataframe,'\n',data)

    for position in sequence_dataframe:
        print(sequence_dataframe[position].value_counts(normalize=True))
        
    return

def oberserved_heterozygosity():
        
    return

def inbreeding_coefficient(sequence_dataframe,data):
    print("calculating inbreeding coefficient")

    print(data)
    print(sequence_dataframe)
    individuals = data['IDV'].unique()
    heteroz_dict = {}
    
    for indv in individuals:
        col_to_concatenate = data['IDV'] == indv
        print(sequence_dataframe[col_to_concatenate]) 
        
        individual_dataframe = sequence_dataframe[col_to_concatenate]    
        
        heteroz_dict[indv] = []
        
        for allele in individual_dataframe:
            print(individual_dataframe[allele][0], individual_dataframe[allele][1])
            allele_combination = ''.join(individual_dataframe[allele].astype('str'))
            heteroz_dict[indv].append(allele_combination)

    diploid_df = pd.DataFrame.from_dict(heteroz_dict)        
    print(diploid_df)
    
    for position in sequence_dataframe:
        frequencies = sequence_dataframe[position].value_counts(normalize=True).tolist()

        H_exp = 1
        for site_freq in frequencies: 
            H_exp *= site_freq
        print(frequencies,H_exp) 
        
    return

def wright_fst(sequence_dataframe,data): 
    print("caluculating wright_fst")   

    return





