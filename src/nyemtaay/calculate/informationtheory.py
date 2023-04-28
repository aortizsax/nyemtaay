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

# Libraries
import numpy as np
import matplotlib.pyplot as plt
import dendropy
import pandas as pd
from nyemtaay.mathlib import mathfn
import os
import math


def sequences_to_gamete_prob(sequence_dataframe, data, identifier,geo):
    print(sequence_dataframe)
    print(data)
    print(identifier)
    
    #convert to gamete probablity dictionary per deme
    
    
    # identififer for subpopulaitons
    ID = identifier

    # initialize denom and numerator
    H_s_sub = {}
    H_i_sub = {}

    # identify subpopulations
    subpopulations = data[ID].unique()

    # H_s
    gamete_dict = {}
    population_dict = {}
    for subpopulation in subpopulations:
        H_s_sub[subpopulation] = []  # positional array
        isSUBPOPULATION = data[ID] == subpopulation
        population_dict[subpopulation] = isSUBPOPULATION.sum()
        subpopulation_sequences = sequence_dataframe[isSUBPOPULATION]
        gamete_dict[subpopulation] = pd.Series(subpopulation_sequences.values.tolist()).map(lambda x: ''.join(map(str,x)))

    gamete_df = pd.DataFrame.from_dict(gamete_dict)
    gamete_counts = gamete_df.apply(pd.Series.value_counts).fillna(0)
    #gamete_sums
    gamete_probabilty = gamete_counts.div(isSUBPOPULATION.sum(),axis=1)
    print(gamete_probabilty)
    
    return (gamete_probabilty, population_dict)
    
    
    
def demes_shannon_entropy(gamete_probabilty):  
    print("shannon_entropy",gamete_probabilty.apply(shannon_entropy,axis=0))
    return gamete_probabilty.apply(shannon_entropy,axis=0)    
    

def shannon_entropy(row_p):
    
    H = -1
    H_sum = 0

    for p_i in row_p:
        H_sum = p_i * np.log(p_i)

    H *= H_sum
    
    return H
    
def demes_jsd(gamete_probabilty, population_dict):
    """
    JSD(P||Q) = ...
    """
    jsd_dict = {}
    D_dict = {}
    I_dict = {}
    for (deme_p, row_p) in gamete_probabilty.items():
        n_p = population_dict[deme_p]
        for (deme_q, row_q) in gamete_probabilty.items():
            n_q = population_dict[deme_q]
            if deme_p != deme_q:
                comparison = deme_p + '->' + deme_q
                print(comparison)
                (w_p, w_q) = statisical_weights(n_p, n_q)
                jsd_dict[comparison] = jsd(row_p, row_q, w_p, w_q)
                D_dict[comparison] = jsd_normalized_to_max(row_p, row_q, w_p, w_q)
                
                X = union(row_p, row_q)
                J = intersection(row_p, row_q)
                print("X",X,"J",J)
                D = disjoint(X,J)
                print("D",D)
                I = information_flow_directionality(row_p,row_q,X,J)
                I_dict[comparison] = I 
                
                
    print("JSD",jsd_dict)
    print("D",D_dict)
    return None

def jsd(row_p, row_q, w_p, w_q):
    
    jsdiv = 0
    
    jsdiv += shannon_entropy(w_p*row_p + w_q*row_q)
    
    jsdiv -= w_p * shannon_entropy(row_p)
    jsdiv -= w_p * shannon_entropy(row_q)
    
    return jsdiv

def union(row_p,row_q):

    states_p = np.where(row_p != 0)
    states_q = np.where(row_q != 0)
    X = np.union1d(states_p,states_q)
    
    return X

def intersection(row_p, row_q): #joint dominium

    states_p = np.where(row_p != 0)
    states_q = np.where(row_q != 0)
    J = np.intersect1d(states_p,states_q)
    
    return J

def disjoint(X,J): #disjoint dominium
    return np.setdiff1d(X,J)

def statisical_weights(n_1,n_2):

    n_T = n_1 + n_2
    
    w_1 = (n_1) / (n_T)    
    w_2 = (n_2) / (n_T)    
        
    return (w_1, w_2)

def jsd_normalized_to_max(row_p, row_q, w_p, w_q):

    D = jsd(row_p, row_q, w_p, w_q)
    demonenator = -w_p * np.log(w_p)
    demonenator += -w_q * np.log(w_q)
    D /= demonenator
    return D

def information_flow_directionality(row_p,row_q,X,J):

    (mu_p, mu_q, mu_pq) = index_mu_PQ(row_p,row_q,X,J)

    P_bar_J = row_p[J]/row_p[J].sum()
    Q_bar_J = row_q[J]/row_q[J].sum()
    
    weighted_subset_H_p = shannon_entropy(P_bar_J) / mu_p
    weighted_subset_H_q = shannon_entropy(Q_bar_J) / mu_q

    I = np.sign(weighted_subset_H_q - weighted_subset_H_p)
    print("I",I>0,I)
    return None


def index_mu_PQ(row_p,row_q,X,J):
    
    mu_p = 0
    mu_q = 0
    mu_pq = 0
    
    for i in J:
        mu_p += row_p[i]
        mu_q += row_q[i]
        mu_pq += row_p[i] / row_q[i]

    demon_p = 0
    demon_q = 0
    demon_pq = 0
    
    for i in X:
        demon_p += row_p[i]
        demon_q += row_q[i]
        demon_pq += row_p[i] / row_q[i]

    mu_p /= demon_p
    mu_q /= demon_q
    mu_pq /= demon_pq
    
    return (mu_p, mu_q, mu_pq)

