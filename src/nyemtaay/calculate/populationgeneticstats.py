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
from nyemtaay.mathlib import mathfn

def segregating_sites(sequences,data):
    print("identifying segratating sites")
    index = data.index
    length = len(sequences[0])
    sequence_dataframe = pd.DataFrame(sequences, 
                                      columns = range(1,length+1),
                                      index = index
                                      )
    
    seg_sites = []
    for position in sequence_dataframe:
        frequncies = sequence_dataframe[position].value_counts(normalize=True)
        if len(frequncies) > 1:
            seg_sites.append(position)
            
    print(seg_sites)
    return seg_sites

def number_segregating_sites(sequences,data):
    print("calculating number segratating sites")
    seg_sites = segregating_sites(sequences,data)
    print("number segrating sites:",len(seg_sites))
    return len(seg_sites)
    
def nucleotide_diversity(sequences,data):
    print("calculating nucleotide diversity(pi)")
    index = data.index
    length = len(sequences[0])
    sequence_dataframe = pd.DataFrame(sequences, 
                                      columns = range(1,length+1),
                                      index = index
                                      )
    

    number_sequences = len(sequence_dataframe)
    number_comparisons = mathfn.nchoose2(number_sequences)
    print(number_sequences)
    

def frequency(sequences,data):
    print("calculating frquencies")
    print(sequences,'\n',data)
    
    index = data.index
    length = len(sequences[0])
    sequence_dataframe = pd.DataFrame(sequences, 
                                      columns = range(1,length+1),
                                      index = index
                                      )
    
    print(sequence_dataframe)

    
    for position in sequence_dataframe:
        print(sequence_dataframe[position].value_counts(normalize=True))
        
    
    return


def wright_fst(sequences,data): 
    print("caluculating wright_fst")   
    index = data.index
    length = len(sequences[0])
    sequence_dataframe = pd.DataFrame(sequences, 
                                      columns = range(1,length+1),
                                      index = index
                                      )
    
    #for position in sequence_dataframe:
     #   print(sequence_dataframe[position].value_counts(normalize=True))
        
    
    return























