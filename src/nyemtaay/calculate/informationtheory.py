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


# Libraries
import numpy as np
import matplotlib.pyplot as plt
import dendropy
import pandas as pd
from nyemtaay.mathlib import mathfn
import os
import math


def sequences_to_gamete_prob(sequence_dataframe, data, ID,geo):
    
    #convert to gamete probablity dictionary per deme

    # identify subpopulations
    subpopulations = data[ID].unique()

    gamete_dict = {}
    population_dict = {}
    for subpopulation in subpopulations:
        isSUBPOPULATION = data[ID] == subpopulation
        population_dict[subpopulation] = isSUBPOPULATION.sum()
        subpopulation_sequences = sequence_dataframe[isSUBPOPULATION]
        gamete_dict[subpopulation] = pd.Series(subpopulation_sequences.values.tolist()).map(lambda x: ''.join(map(str,x)))

    gamete_df = pd.DataFrame.from_dict(gamete_dict)
    gamete_df#randomized combinations of this!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    gamete_counts = gamete_df.apply(pd.Series.value_counts).fillna(0)
    #gamete_sums
    gamete_probabilty = gamete_counts.div(isSUBPOPULATION.sum(),axis=1)
    
    return (gamete_probabilty, population_dict)#add random gamete sequences??
    

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
    I_R_dict_list = {}
    I_R_dict_mu_std = {}
    for (deme_p, row_p) in gamete_probabilty.items():
        n_p = population_dict[deme_p]
        for (deme_q, row_q) in gamete_probabilty.items():
            n_q = population_dict[deme_q]
            if deme_p != deme_q:
                comparison = deme_p + '->' + deme_q
                (w_p, w_q) = statisical_weights(n_p, n_q)
                jsd_dict[comparison] = jsd(row_p, row_q, w_p, w_q)
                
        
                
    print("JSD",jsd_dict)

    return jsd_dict
    
    
def jsd(row_p, row_q, w_p, w_q):
    """
    JSD(P||Q) = ...
    """
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

def demes_norm_jsd(gamete_probabilty, population_dict):
    """
    JSD(P||Q) = ...
    """
    jsd_dict = {}
    D_dict = {}
    I_dict = {}
    I_R_dict_list = {}
    I_R_dict_mu_std = {}
    for (deme_p, row_p) in gamete_probabilty.items():
        n_p = population_dict[deme_p]
        for (deme_q, row_q) in gamete_probabilty.items():
            n_q = population_dict[deme_q]
            if deme_p != deme_q:
                comparison = deme_p + '->' + deme_q
                (w_p, w_q) = statisical_weights(n_p, n_q)
                jsd_dict[comparison] = jsd(row_p, row_q, w_p, w_q)
                D_dict[comparison] = jsd_normalized_to_max(row_p, row_q, w_p, w_q)
                
                X = union(row_p, row_q)
                J = intersection(row_p, row_q)
                D = disjoint(X,J)
                I = information_flow_directionality(row_p,row_q,X,J)
                I_dict[comparison] = I 
                
                #shuffle population sequences 1000 times 
                #get norm mean and std
                
                
    print("D",D_dict)

    return D_dict


def jsd_normalized_to_max(row_p, row_q, w_p, w_q):

    D = jsd(row_p, row_q, w_p, w_q)
    demonenator = -w_p * np.log(w_p)
    demonenator += -w_q * np.log(w_q)
    D /= demonenator
    return D
    
    
def demes_info_flow_direction(gamete_probabilty, population_dict):
    """
    I(P||Q) = ...
    """
    I_dict = {}
    I_R_dict_list = {}
    I_R_dict_mu_std = {}
    for (deme_p, row_p) in gamete_probabilty.items():
        n_p = population_dict[deme_p]
        for (deme_q, row_q) in gamete_probabilty.items():
            n_q = population_dict[deme_q]
            if deme_p != deme_q:
                comparison = deme_p + '->' + deme_q
                (w_p, w_q) = statisical_weights(n_p, n_q)
                
                
                
                X = union(row_p, row_q)
                J = intersection(row_p, row_q)
                D = disjoint(X,J)
                I = information_flow_directionality(row_p,row_q,X,J)
                I_dict[comparison] = I 
                
                #shuffle population sequences 1000 times 
                #get norm mean and std
                

    print("I",I_dict)
    return I_dict

def information_flow_directionality(row_p,row_q,X,J):
    """
    """
    (mu_p, mu_q, mu_pq) = index_mu_PQ(row_p,row_q,X,J)

    P_bar_J = row_p[J]/row_p[J].sum()
    Q_bar_J = row_q[J]/row_q[J].sum()
    
    weighted_subset_H_p = shannon_entropy(P_bar_J) / mu_p
    weighted_subset_H_q = shannon_entropy(Q_bar_J) / mu_q

    I = np.sign(weighted_subset_H_q - weighted_subset_H_p)
    return I


def index_mu_PQ(row_p,row_q,X,J):
    """
    """
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



def sequences_to_random_deme_combinations(sequence_dataframe, data, identifier,geo):
    
    #convert to gamete probablity dictionary per deme
    rng = np.random.default_rng(12345)
    
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
    gamete_df#randomized combinations of this!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    gamete_random_combinations_dict = {}
    for (deme_p,gamete_sequences) in gamete_df.items():

        for (deme_q, row_q) in gamete_df.items():
            if deme_p != deme_q:
                comparison = deme_p + '->' + deme_q
                gamete_random_combinations_dict[comparison] = []

                for i in range(100):
                    comparison_sequences = gamete_df[[deme_p,deme_q]].copy(deep=True)
                    #.sample(frac=1,random_state=1)
                    
                    
                    dn=[]
                    for deme in [deme_p,deme_q]:
                        list_of_random_sequences = []
                        
                        col = rng.integers(low=0, high=2, size=population_dict[deme_p])
                        row = rng.integers(low=0, high=population_dict[deme_p], size=population_dict[deme_p])
                        
                        for element_index, random_index in enumerate(col):
                            random_jndex = row[element_index]
                            list_of_random_sequences.append(comparison_sequences.iloc[:,random_index][random_jndex])

                        df1 = pd.DataFrame({'sequence': range(population_dict[deme]),
                                     deme: list_of_random_sequences,
                                    }).set_index('sequence')
                        dn.append(df1)
                            
                    comparison_random_sequences = pd.concat(dn, axis=1)
                    
                    rand_gamete_counts = comparison_random_sequences.apply(pd.Series.value_counts).fillna(0)
                    #gamete_sums
                    rand_gamete_probabilty = rand_gamete_counts.div(isSUBPOPULATION.sum(),axis=1)                

                    gamete_random_combinations_dict[comparison].append(rand_gamete_probabilty)
    
    
    return (gamete_random_combinations_dict, population_dict)#add random gamete sequences??

def randomized_information_flow_directionality(gamete_random_combinations_dict,I_dict):
    I_R_array = []
    
    IR_dict ={}
    R_dict ={}
    
    for comparison, list_of_gamete_probabilites in gamete_random_combinations_dict.items(): 
        deme_p = comparison[0]
        deme_q = comparison[-1]
        I_pq = I_dict[comparison]
        
        
        for i in range(100):
            gamete_probabilty = list_of_gamete_probabilites[i] 
            row_p = gamete_probabilty[deme_p]
            row_q = gamete_probabilty[deme_q]
            
            X = union(row_p, row_q) 
            J = intersection(row_p, row_q) 
            D = disjoint(X,J) 
            I = information_flow_directionality(row_p,row_q,X,J) 
            I_R_array.append(I)


        std_I_R = np.std(I_R_array, dtype=np.float64)
        mu_I_R = np.mean(I_R_array, dtype=np.float64)
        IR_dict[comparison] = (mu_I_R, std_I_R)
        
        R = abs(I_pq - mu_I_R)
        R /= std_I_R
        
        R_dict[comparison] = R
        
    print("IR ",IR_dict)
    print("R ",R_dict)    
    
        
    return IR_dict

