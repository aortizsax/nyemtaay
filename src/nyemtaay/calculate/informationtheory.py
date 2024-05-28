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
#get directionality of I on to network

# Libraries
import numpy as np
import matplotlib.pyplot as plt
import dendropy
import pandas as pd
from nyemtaay.mathlib import mathfn
from nyemtaay.cluster import informationtheoryclustering
import os
import math


def sequences_to_gamete_prob(sequence_dataframe, data, ID,geo):
    #print(sequence_dataframe)
    #convert to gamete probablity dictionary per deme

    # identify subpopulations
    subpopulations = data[ID].unique()

    gamete_dict = {}
    allele_dict = {}
    population_dict = {}
    pos_allele_probablity_deme_dict = {}
    for subpopulation in subpopulations:
        isSUBPOPULATION = data[ID] == subpopulation
        population_dict[subpopulation] = isSUBPOPULATION.sum()
        subpopulation_sequences = sequence_dataframe[isSUBPOPULATION]
        allele_dict[subpopulation] = subpopulation_sequences
        gamete_dict[subpopulation] = pd.Series(subpopulation_sequences.values.tolist()).map(lambda x: ''.join(map(str,x)))

    gamete_df = pd.DataFrame.from_dict(gamete_dict)
    gamete_counts = gamete_df.apply(pd.Series.value_counts).fillna(0)
    #print(gamete_counts)
    gamete_probabilty = {}
    for subpopulation in subpopulations:
        isSUBPOPULATION = data[ID] == subpopulation

        gamete_probabilty[subpopulation] = gamete_counts[subpopulation] / isSUBPOPULATION.sum()

    gamete_probabilty = pd.DataFrame.from_dict(gamete_probabilty)
    #print(gamete_probabilty)
    for key in allele_dict.keys():
        allele_df = allele_dict[key]
        allele_counts = allele_df.apply(pd.Series.value_counts).fillna(0)
        allele_probabilty = allele_counts.div(isSUBPOPULATION.sum(),axis=1)
        pos_allele_probablity_deme_dict[key] = allele_probabilty

    return (gamete_probabilty, population_dict,pos_allele_probablity_deme_dict) #add random gamete sequences??



def demes_shannon_entropy(gamete_probabilty,geo):
    shannon_entropy_by_deme = gamete_probabilty.apply(shannon_entropy,axis=0)
#    print( gamete_probabilty)
#    print(shannon_entropy_by_deme)
    weighted_bool = False
    (networkx_dictionary,node_color_array) = information_theory_dict_to_pd_df(shannon_entropy_by_deme,
                                                                    weighted_bool,
                                                                    geo,
                                                                    'node')
    # plot_node_metric(networkx_dictionary,node_color_array)

    # print("shannon_entropy by node",shannon_entropy_by_deme)
    return (shannon_entropy_by_deme , False)

def demes_allele_shannon_entropy(pos_allele_probablity_deme_dict):
    for deme in pos_allele_probablity_deme_dict.keys():
        allele_frequencies = pos_allele_probablity_deme_dict[deme]
        print("shannon_entropy",deme,allele_frequencies.apply(shannon_entropy,axis=0))

    return None

def shannon_entropy(row_p):
    H = 0

    for p_i in row_p:
        if p_i > 0:
            H -= p_i * np.log(p_i) #dont forget +=

    return H

def demes_jsd(gamete_probabilty, population_dict,geo):
    """
    JSD(P||Q) = ...
    """
    jsd_dict = {}
    for (deme_p, row_p) in gamete_probabilty.items():
        n_p = population_dict[deme_p]
        for (deme_q, row_q) in gamete_probabilty.items():
            n_q = population_dict[deme_q]
            if deme_p != deme_q:
                comparison = deme_p + '->' + deme_q
                (w_p, w_q) = statisical_weights(n_p, n_q)
                jsd_dict[comparison] = jsd(row_p, row_q, w_p, w_q)



    print("JSD unidirectional edge",jsd_dict)

    weighted_bool = False
    (networkx_dictionary,node_color_array) = information_theory_dict_to_pd_df(jsd_dict,
                                                                    weighted_bool,
                                                                    geo,
                                                                    'edge')
    percolation_network = plot_unidirectional_metric(networkx_dictionary,node_color_array)
    return (jsd_dict, percolation_network)

def demes_allele_jsd(pos_allele_probablity_deme_dict, population_dict):
    """
    """
    jsd_dict = {}
    for (deme_p, row_pp) in pos_allele_probablity_deme_dict.items():

        n_p = population_dict[deme_p]
        for (deme_q, row_qq) in pos_allele_probablity_deme_dict.items():
            n_q = population_dict[deme_q]
            if deme_p != deme_q:
                for i, row_p in row_pp.items():
                    row_q = row_qq[i]
                    comparison = deme_p + '->' + deme_q + '_' + str(i)
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

    jsdiv -= w_q * shannon_entropy(row_q)

    return jsdiv

def union(row_p,row_q):

    states_p = np.where(row_p > 0)
    states_q = np.where(row_q > 0)
    X = np.union1d(states_p,states_q)

    return X

def intersection(row_p, row_q): #joint dominium

    states_p = np.where(row_p > 0)
    states_q = np.where(row_q > 0)
    J = np.intersect1d(states_p,states_q)

    return J

def disjoint(X,J): #disjoint dominium
    return np.setdiff1d(X,J)

def statisical_weights(n_1,n_2):

    n_T = n_1 + n_2

    w_1 = (n_1) / (n_T)
    w_2 = (n_2) / (n_T)

    return (w_1, w_2)

def demes_norm_jsd(gamete_probabilty, metadata, population_dict,geo):
    """
    JSD(P||Q) = ...
    """
    D_dict = {}
    for i , (deme_p, row_p) in enumerate(gamete_probabilty.items()):
        n_p = population_dict[deme_p]
        for j, (deme_q, row_q) in enumerate(gamete_probabilty.items()):
            n_q = population_dict[deme_q]
            if j > i: #if deme_p != deme_q:
                comparison = deme_p + '->' + deme_q
                (w_p, w_q) = statisical_weights(n_p, n_q)

                D_dict[comparison] = jsd_normalized_to_max(row_p, row_q, w_p, w_q)

    weighted_bool = False
    print("D",D_dict)
    (networkx_dictionary,node_color_array) = information_theory_dict_to_pd_df(D_dict,
                                                                    weighted_bool,
                                                                    geo,
                                                                    'edge')
    (percolation_network, H, dist_matrix) = plot_unidirectional_metric(networkx_dictionary,node_color_array)


    informationtheoryclustering.dbscan_imp(dist_matrix)
    # informationtheoryclustering.louvian_clustering(H,metadata)
    return (D_dict, percolation_network)

def demes_allele_norm_jsd(pos_allele_probablity_deme_dict, population_dict): ####################################
    """
    D(P||Q) = JSD / max
    """
    D_dict = {}
    for (deme_p, row_pp) in pos_allele_probablity_deme_dict.items():
        n_p = population_dict[deme_p]
        for (deme_q, row_qq) in pos_allele_probablity_deme_dict.items():
            n_q = population_dict[deme_q]
            if deme_p != deme_q:
                for i, row_p in row_pp.items():
                    row_q = row_qq[i]
                    comparison = deme_p + '->' + deme_q + '_' + str(i)
                    (w_p, w_q) = statisical_weights(n_p, n_q)

                    D_dict[comparison] = jsd_normalized_to_max(row_p, row_q, w_p, w_q)



    print("D",D_dict)

    return D_dict


def jsd_normalized_to_max(row_p, row_q, w_p, w_q):

    D = jsd(row_p, row_q, w_p, w_q)
    demonenator = -w_p * np.log(w_p)
    demonenator -= w_q * np.log(w_q)
    D /= demonenator
    return D


def demes_info_flow_direction(gamete_probabilty, population_dict,perc_network,geo):
    """
    I(P||Q) = ...
    """
    I_dict = {}
    w_dict = {}
    I_R_dict_list = {}
    I_R_dict_mu_std = {}
    for index, row in perc_network.iterrows():
        deme_p = row['Source']
        deme_q = row['Target']
        row_p = gamete_probabilty[deme_p]
        row_q = gamete_probabilty[deme_q]
        n_p = population_dict[deme_p]
        n_q = population_dict[deme_q]
        comparison = deme_p + '->' + deme_q


        X = union(row_p, row_q)
        # print('X',comparison,X)
        J = intersection(row_p, row_q)
        # print('J',comparison,J)#w_p,w_q)
        D = disjoint(X,J)
        # print('D',comparison,D)
       # print(row_p,row_q)
        (mu_p, mu_q, mu_pq) = index_mu_PQ(row_p,row_q,X,J)
        #print(comparison,mu_p,mu_q)
        I = information_flow_directionality(row_p,row_q,X,J)
        I_dict[comparison] = I
        w_dict[comparison] = (mu_p, mu_q)
#    for (deme_p, row_p) in gamete_probabilty.items():
#        n_p = population_dict[deme_p]
#        for (deme_q, row_q) in gamete_probabilty.items():
#            n_q = population_dict[deme_q]
#            if deme_p != deme_q:
#                comparison = deme_p + '->' + deme_q
#                (w_p, w_q) = statisical_weights(n_p, n_q)
#
#
#
#                X = union(row_p, row_q)
#                J = intersection(row_p, row_q)
#                D = disjoint(X,J)
#                I = information_flow_directionality(row_p,row_q,X,J)
#                I_dict[comparison] = I

                #shuffle population sequences 1000 times
                #get norm mean and std


    # print("I",I_dict)
    # print('w',w_dict)
    return (I_dict , w_dict)

def demes_allele_info_flow_direction(pos_allele_probablity_deme_dict, population_dict,perc_network):
    """
    I(P||Q) = ...
    """
    print(perc_network)
    exit()
    I_dict = {}
    w_dict = {}
    I_R_dict_list = {}
    I_R_dict_mu_std = {}

    for (deme_p, row_pp) in pos_allele_probablity_deme_dict.items():
        n_p = population_dict[deme_p]
        for (deme_q, row_qq) in pos_allele_probablity_deme_dict.items():
            n_q = population_dict[deme_q]
            if deme_p != deme_q:
                for i, row_p in row_pp.items():
                    row_q = row_qq[i]
                    comparison = deme_p + '->' + deme_q + '_' + str(i)
                    (w_p, w_q) = statisical_weights(n_p, n_q)



                    X = union(row_p, row_q)
                    J = intersection(row_p, row_q)

                    D = disjoint(X,J)
                    print(X,J,D)
                    I = information_flow_directionality(row_p,row_q,X,J)
                    I_dict[comparison] = I
                    w_dict[comparison] = (w_p, w_q)
                #shuffle population sequences 1000 times
                #get norm mean and std


    print("I",I_dict)
    print('w',w_dict)
    return (I_dict , w_dict)

def information_flow_directionality(row_p,row_q,X,J): #pass index mu to this to call less
    """
    """
    (mu_p, mu_q, mu_pq) = index_mu_PQ(row_p,row_q,X,J)


    P_bar_J = row_p[J]#/row_p[J].sum()
    Q_bar_J = row_q[J]#/row_q[J].sum()

    weighted_subset_H_p = shannon_entropy(P_bar_J) / mu_p
    weighted_subset_H_q = shannon_entropy(Q_bar_J) / mu_q

    I = np.sign(weighted_subset_H_q - weighted_subset_H_p)
    return I


def index_mu_PQ(row_p,row_q,X,J):
    """
    """
#    print(X,J)
    mu_p = 0.0
    mu_q = 0.0
    mu_pq = 0.0

    for i in J:
        mu_p += row_p[i]
        mu_q += row_q[i]
        #mu_pq += row_p[i] / row_q[i]

#    demon_p = 0.0
#    demon_q = 0.0
#    demon_pq = 0.0
#
#    for i in X:
#        demon_p += row_p[i]
#        demon_q += row_q[i]
#        #demon_pq += row_p[i] / row_q[i]

#    mu_p /= demon_p
#    mu_q /= demon_q
    #mu_pq /= demon_pq
#    print(mu_p, mu_q, mu_pq)
    return (mu_p, mu_q, mu_pq)



def sequences_to_random_deme_combinations(sequence_dataframe, data, identifier,perc_network,geo):

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
    max_pop_size = max(list(population_dict.values()))
    gamete_df = pd.DataFrame.from_dict(gamete_dict)
    gamete_df#randomized combinations of this!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    gamete_random_combinations_dict = {}
    for index, row in perc_network.iterrows():
        deme_p = row['Source']
        deme_q = row['Target']
        gamete_sequences = gamete_df[deme_p]
        row_q = gamete_df[deme_q]
        n_p = population_dict[deme_p]
        n_q = population_dict[deme_q]
        comparison = deme_p + '->' + deme_q
        gamete_random_combinations_dict[comparison] = []
        for i in range(1000):

            comparison_sequences = gamete_df[[deme_p,deme_q]].copy(deep=True)
            #.sample(frac=1,random_state=1)


            dn=[]
            for deme in [deme_p,deme_q]:
                list_of_random_sequences = []

                col = rng.integers(low=0, high=2, size=max_pop_size)
                row = rng.integers(low=0, high=max_pop_size, size=max_pop_size)

                for element_index, random_index in enumerate(col):

                    random_jndex = row[element_index]%population_dict[deme]
                    list_of_random_sequences.append(comparison_sequences.iloc[:,random_index][random_jndex])

                df1 = pd.DataFrame({'sequence': range(max_pop_size),
                             deme: list_of_random_sequences,
                            }).set_index('sequence')
                dn.append(df1)

            comparison_random_sequences = pd.concat(dn, axis=1)

            rand_gamete_counts = comparison_random_sequences.apply(pd.Series.value_counts).fillna(0)
            #gamete_sums
            rand_gamete_probabilty = rand_gamete_counts.div(isSUBPOPULATION.sum(),axis=1)

            gamete_random_combinations_dict[comparison].append(rand_gamete_probabilty)


    return (gamete_random_combinations_dict, population_dict)#add random gamete sequences??

def randomized_information_flow_directionality(gamete_random_combinations_dict,I_dict, geo, plot_output_path):
    I_R_array = []

    IR_dict ={}
    R_dict ={}

    for comparison, list_of_gamete_probabilites in gamete_random_combinations_dict.items():

        deme_p = comparison.split('->')[0]
        deme_q = comparison.split('->')[-1]
        I_pq = I_dict[comparison]

        for i in range(1000):

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

    weighted_bool = True
    (networkx_dictionary,node_color_array) = directionality_dicts_to_pd_df(R_dict,
                                                                    I_dict,
                                                                    weighted_bool,
                                                                    geo,
                                                                    'edge')
    G = plot_bidirectional_metric(networkx_dictionary, node_color_array, plot_output_path)

    #informationtheoryclustering.louvian_clustering(G)

    return (R_dict, weighted_bool)

def directionality_dicts_to_pd_df(information_theory_dictionary,direction_dict,weighted_bool,geometry,metric_char):
    """
    """
    #weighted bool == True
    threshold = -1110.023
    if weighted_bool:
        networkx_dictionary = {"Source":[],"Target":[],"Type":[],"weight":[],"mass":[]}
        node_color_array = []
        for comparison, edge_weight in information_theory_dictionary.items():
            #print(comparison)
            if direction_dict[comparison] == -1:
                deme_source = comparison.split('->')[-1]
                deme_target = comparison.split('->')[0]

            elif direction_dict[comparison] >= 0:
                deme_source = comparison.split('->')[0]
                deme_target = comparison.split('->')[-1]
            node_color_array.append(1)
           # for i, deme_target in enumerate(fxt_dictionary.keys()):
            #    if deme_source != deme_target:
            networkx_dictionary["Source"].append(deme_source)
            networkx_dictionary["Target"].append(deme_target)
            networkx_dictionary["Type"].append("Directed")
            if str(edge_weight) == 'nan':
                edge_weight = 0
            networkx_dictionary["weight"].append(str(round(edge_weight,3)))
            networkx_dictionary["mass"].append(1)

    else:
        if metric_char == "edge":
            networkx_dictionary = {"Source":[],"Target":[],"Type":[],"weight":[],"mass":[]}
            node_color_array = []
            for comparison, edge_weight in information_theory_dictionary.items():
                if edge_weight > threshold:
                    #print(comparison)
                    deme_source = comparison.split('->')[0]
                    deme_target = comparison.split('->')[-1]
                    node_color_array.append(1)
                   # for i, deme_target in enumerate(fxt_dictionary.keys()):
                    #    if deme_source != deme_target:
                    networkx_dictionary["Source"].append(deme_source)
                    networkx_dictionary["Target"].append(deme_target)
                    networkx_dictionary["Type"].append("Unidirectional")
                    networkx_dictionary["weight"].append(round(edge_weight,3))
                    networkx_dictionary["mass"].append(1)

        elif metric_char == "node":
            networkx_dictionary = {"Source":[],"Target":[],"Type":[],"weight":[],"mass":[]}
            node_color_array = []
            for i, (deme_source, mass) in enumerate(information_theory_dictionary.items()):
                node_color_array.append(mass)
                if geometry == 'complete graph':
                    for j, deme_target in enumerate(information_theory_dictionary.keys()):
                        if deme_source != deme_target:
                            networkx_dictionary["Source"].append(deme_source)
                            networkx_dictionary["Target"].append(deme_target)
                            networkx_dictionary["Type"].append("Undirected")
                            networkx_dictionary["weight"].append(1)
                            networkx_dictionary["mass"].append(mass)

                elif geometry == 'chain graph':

                    if i ==0 :
                        networkx_dictionary["Source"].append(deme_source)
                        networkx_dictionary["Target"].append(list(information_theory_dictionary.keys())[i+1])
                        networkx_dictionary["Type"].append("Undirected")
                        networkx_dictionary["weight"].append(1)
                        networkx_dictionary["mass"].append(mass)


                    elif i == len(list(information_theory_dictionary.keys())) -1:
                        pass
                        networkx_dictionary["Source"].append(deme_source)
                        networkx_dictionary["Target"].append(list(information_theory_dictionary.keys())[i-1])
                        networkx_dictionary["Type"].append("Undirected")
                        networkx_dictionary["weight"].append(1)
                        networkx_dictionary["mass"].append(mass)

                    else :
                        networkx_dictionary["Source"].append(deme_source)
                        networkx_dictionary["Target"].append(list(information_theory_dictionary.keys())[i+1])
                        networkx_dictionary["Type"].append("Undirected")
                        networkx_dictionary["weight"].append(1)
                        networkx_dictionary["mass"].append(mass)

#        print(networkx_dictionary)
    return (networkx_dictionary,node_color_array)



def information_theory_dict_to_pd_df(information_theory_dictionary,weighted_bool,geometry,metric_char):
    """
    """
    #weighted bool == True
    threshold = -1110.023
    if weighted_bool:
        networkx_dictionary = {"Source":[],"Target":[],"Type":[],"weight":[],"mass":[]}
        node_color_array = []
        for comparison, edge_weight in information_theory_dictionary.items():
            #print(comparison)
            deme_source = comparison.split('->')[0]
            deme_target = comparison.split('->')[-1]
            node_color_array.append(1)
           # for i, deme_target in enumerate(fxt_dictionary.keys()):
            #    if deme_source != deme_target:
            networkx_dictionary["Source"].append(deme_source)
            networkx_dictionary["Target"].append(deme_target)
            networkx_dictionary["Type"].append("Directed")
            networkx_dictionary["weight"].append(str(round(edge_weight,3)))
            networkx_dictionary["mass"].append(1)

    else:
        if metric_char == "edge":
            networkx_dictionary = {"Source":[],"Target":[],"Type":[],"weight":[],"mass":[]}
            node_color_array = []
            for comparison, edge_weight in information_theory_dictionary.items():
                if edge_weight > threshold:
                    #print(comparison)
                    deme_source = comparison.split('->')[0]
                    deme_target = comparison.split('->')[-1]
                    node_color_array.append(1)
                   # for i, deme_target in enumerate(fxt_dictionary.keys()):
                    #    if deme_source != deme_target:
                    networkx_dictionary["Source"].append(deme_source)
                    networkx_dictionary["Target"].append(deme_target)
                    networkx_dictionary["Type"].append("Unidirectional")
                    networkx_dictionary["weight"].append(round(edge_weight,3))
                    networkx_dictionary["mass"].append(1)

        elif metric_char == "node":
            networkx_dictionary = {"Source":[],"Target":[],"Type":[],"weight":[],"mass":[]}
            node_color_array = []
            for i, (deme_source, mass) in enumerate(information_theory_dictionary.items()):
                node_color_array.append(mass)
                if geometry == 'complete graph':
                    for j, deme_target in enumerate(information_theory_dictionary.keys()):
                        #if j > i:
                        if deme_source != deme_target:
                            networkx_dictionary["Source"].append(deme_source)
                            networkx_dictionary["Target"].append(deme_target)
                            networkx_dictionary["Type"].append("Undirected")
                            networkx_dictionary["weight"].append(1)
                            networkx_dictionary["mass"].append(mass)

                elif geometry == 'chain graph':

                    if i ==0 :
                        networkx_dictionary["Source"].append(deme_source)
                        networkx_dictionary["Target"].append(list(information_theory_dictionary.keys())[i+1])
                        networkx_dictionary["Type"].append("Undirected")
                        networkx_dictionary["weight"].append(1)
                        networkx_dictionary["mass"].append(mass)


                    elif i == len(list(information_theory_dictionary.keys())) -1:
                        pass
                        networkx_dictionary["Source"].append(deme_source)
                        networkx_dictionary["Target"].append(list(information_theory_dictionary.keys())[i-1])
                        networkx_dictionary["Type"].append("Undirected")
                        networkx_dictionary["weight"].append(1)
                        networkx_dictionary["mass"].append(mass)

                    else :
                        networkx_dictionary["Source"].append(deme_source)
                        networkx_dictionary["Target"].append(list(information_theory_dictionary.keys())[i+1])
                        networkx_dictionary["Type"].append("Undirected")
                        networkx_dictionary["weight"].append(1)
                        networkx_dictionary["mass"].append(mass)

#        print(networkx_dictionary)
    return (networkx_dictionary,node_color_array)




#################################################################################

def plot_node_metric(networkx_format_dictionary,node_color_array):
    """
    """
    import networkx as nx
    import matplotlib.pyplot as plt

    networkx_df = pd.DataFrame.from_dict(networkx_format_dictionary)
##    print("bynode")
    G = nx.from_pandas_edgelist(networkx_df,
                                source="Source",
                                target="Target",
                                edge_attr="weight")

    pos = nx.circular_layout(G)



    fig, ax = plt.subplots()
    nx.draw_networkx_nodes(G,
                            pos,
                            node_color=node_color_array, #set floor and ceiling
                            node_size=800, # popsize in the future
                            cmap=plt.cm.Blues,
                            ax=ax, vmin=0, vmax=1)
    nx.draw_networkx_labels(G, pos, ax=ax)

    straight_edges = list(set(G.edges()))
    nx.draw_networkx_edges(G, pos, ax=ax, edgelist=straight_edges)

    plt.show()

    return None

def flatten_dict(dict_of_dicts):
    # Flatten Nested Dictionary to Matrix using pandas DataFrame
    df = pd.DataFrame(dict_of_dicts).T.fillna(0)

    # printing result
    res = df.reset_index().rename(columns={'index': 'Grouped'})
    res = [list(res.columns)] + res.values.tolist()
    return res

def plot_unidirectional_metric(networkx_format_dictionary,node_color_array):

    """
    """
    import networkx as nx
    import matplotlib.pyplot as plt
    networkx_df = pd.DataFrame.from_dict(networkx_format_dictionary)
    #print(networkx_df)
    print("byunidedge")
    G = nx.from_pandas_edgelist(networkx_df,
                                source="Source",
                                target="Target",
                                edge_attr="weight"
                                )
    dist_matrix = dict(nx.all_pairs_dijkstra_path_length(G))
    #print('matrix',np.delete(np.delete(np.array(flatten_dict(dist_matrix)),0,0),0,1))
    dist_matrix = np.delete(np.delete(np.array(flatten_dict(dist_matrix)),0,0),0,1)
    dist_matrix = np.square(dist_matrix.astype(float))

    #G = nx.Graph(networkx_format_dictionary)
    #node_color_array = networkx_format_dictionary['mass']

    pos = nx.circular_layout(G)

    T=nx.minimum_spanning_tree(G)
    # plt.figure()
    # nx.draw_networkx(T, pos=pos, with_labels=True, node_size = 500)
    # plt.show()

    threshold = 0.02

    distmatrix = nx.attr_matrix(T,edge_attr="weight")[0]
    #print("TYPE",G.edges(data=True),type(G.edges(data=True)))
    threshold = find_threshold_bfs(distmatrix)
    threshold = distmatrix.max()
    print('Threshold is {}'.format(threshold))


    networkx_df = pd.DataFrame.from_dict(networkx_format_dictionary)
    networkx_df_thresh = networkx_df.loc[lambda networkx_df: networkx_df['weight'] <= threshold, :]

    #networkx_df_thresh = nx.to_pandas_edgelist(T, nodelist=list(T))  #maybe include ekey
    print("byunidedge")
    G = nx.from_pandas_edgelist(networkx_df_thresh,
                                source="Source",
                                target="Target",
                                edge_attr="weight"
                                )
    H = nx.Graph()
    H.add_nodes_from(sorted(G.nodes(data=True)))
    H.add_edges_from(G.edges(data=True))


    pos = nx.circular_layout(H)

    fig, ax = plt.subplots()
    nx.draw_networkx_nodes(H,
                            pos,
                            ax=ax)
    nx.draw_networkx_labels(H, pos, ax=ax)

    straight_edges = np.array(list(set(H.edges())))
    edge_labels = nx.get_edge_attributes(H, "weight")
    print("TYPE",edge_labels,type(edge_labels))

    # nx.draw_networkx_edges(H, pos, ax=ax, edgelist=straight_edges)

    # plt.show()

    d = nx.to_dict_of_dicts(H)
    #print(d)
    for k,v in d.items():
        for k2,v2 in v.items():
            d[k][k2]['weight'] = v[k2]['weight']**2
    H = nx.Graph(d)

    return (networkx_df_thresh, H, dist_matrix)

def find_threshold_bfs(array):
    #print(type(array))
    first_node = 0
    last_node = len(array) - 1
    probabilities = np.unique(list(array.flat))#flatten
    low = 0
    high = len(probabilities)

    while high - low > 1:
        i = (high + low) // 2
        prob = probabilities[i]
        copied_array = np.array(array)
        copied_array[copied_array < prob] = 0.0
        if bfs(copied_array, first_node, last_node):
            low = i
        else:
            high = i

    return probabilities[low]


def bfs(graph, source, dest):
    """Perform breadth-first search starting at source. If dest is reached,
    return True, otherwise, return False."""
    # Based on http://www.ics.uci.edu/~eppstein/PADS/BFS.py
    # by D. Eppstein, July 2004.
    visited = set([source])
    nodes = np.arange(0, len(graph))
    stack = [(source, nodes[graph[source] > 0])]
    while stack:
        parent, children = stack[0]
        for child in children:
            if child == dest:
                return True
            if child not in visited:
                visited.add(child)
                stack.append((child, nodes[graph[child] > 0]))
        stack.pop(0)
    return False




def plot_bidirectional_metric(networkx_format_dictionary, node_color_array, plot_output_path):
    """
    """
    import networkx as nx
    import matplotlib.pyplot as plt

    networkx_df = pd.DataFrame.from_dict(networkx_format_dictionary)
    networkx_df = networkx_df.fillna(0)

    G = nx.from_pandas_edgelist(networkx_df,
                                source="Source",
                                target="Target",
                                edge_attr="weight",
                                create_using=nx.DiGraph())

    edge_list = nx.to_pandas_edgelist(G)

    pos = nx.circular_layout(G)

    cmap = plt.cm.plasma


    H = nx.DiGraph()
    H.add_nodes_from(sorted(G.nodes(data=True)))
    H.add_weighted_edges_from(G.edges(data=True))

    pos = nx.circular_layout(H)
    nx.draw_networkx(H, pos)

    for edge in H.edges(data='weight'):
        #print(edge,dir(edge))
        weight = float(edge[2]['weight'])
        thickness = weight
        if weight > 1:
            thickness = weight * weight
        #print(thickness)
        nx.draw_networkx_edges(H, pos, edgelist=[edge], width=thickness,label=str(round(weight,2)))
    # plt.show()
    plot.savefig(plot_output_path)
    plot.close()

    return H





