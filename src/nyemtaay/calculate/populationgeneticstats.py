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
# equations from https://www.biorxiv.org/content/10.1101/2020.01.30.927186v1.full.pdf
# and https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8527463/
#
###########################################################################
## Reader Implementation

# Libraries
import numpy as np
import matplotlib.pyplot as plt
import dendropy
import pandas as pd
from nyemtaay.mathlib import mathfn
import os
import math


def segregating_sites(sequence_dataframe, data):
    print("identifying segratating sites")

    seg_sites = []
    for position in sequence_dataframe:
        frequncies = sequence_dataframe[position].value_counts(normalize=True)
        if len(frequncies) > 1:
            seg_sites.append(position)

    print(seg_sites)
    return seg_sites


def number_segregating_sites(sequence_dataframe, data):
    print("calculating number segratating sites")
    seg_sites = segregating_sites(sequence_dataframe, data)
    print("number segrating sites:", len(seg_sites))
    return len(seg_sites)


def number_pairwise_differences(sequence_dataframe):
    count_pairwise_differences = 0

    for i, (index, sequence) in enumerate(sequence_dataframe[:-1].iterrows()):
        # print(i,sequence)
        difference = abs(sequence - sequence_dataframe[i + 1 :])
        # print(difference)
        normalized_diff = difference / difference

        count_pairwise_differences += normalized_diff.sum().sum()
    return count_pairwise_differences


def nucleotide_diversity(sequence_dataframe, data):
    print("calculating nucleotide diversity(pi)")

    number_sequences = len(sequence_dataframe)
    number_comparisons = mathfn.nchoose2(number_sequences)
    length = len(sequence_dataframe.iloc[0])

    nucleotide_div = number_pairwise_differences(sequence_dataframe)
    nucleotide_div /= number_comparisons
    nucleotide_div /= length

    print("nucdiv", nucleotide_div)
    return nucleotide_div


def frequency_spectrum(sequence_dataframe, data):
    print("calculating sfs (add per sequnce length)")
    length = len(sequence_dataframe.iloc[0])
    print("length", length)

    site_frequency_spectrum = [0] * int(length)
    for position in sequence_dataframe:
        print("posisiton", position)
        print(sequence_dataframe[position])
        site_counts = sequence_dataframe[position].value_counts().tolist()
        print(site_counts)
        if len(site_counts) > 1:
            site_frequency_spectrum[site_counts[1]] += 1

    print(site_frequency_spectrum)

    y_pos = np.arange(len(site_frequency_spectrum))

    plot_length = int(len(site_frequency_spectrum))

    plt.bar(y_pos[:plot_length], site_frequency_spectrum[:plot_length])
    plt.show()

    # write to output folder
    path = "./sfs_output"
    # Check whether the specified path exists or not
    isExist = os.path.exists(path)

    # unmask
    os.umask(0)

    # check output folder
    if isExist == False:
        # Create a new directory because it does not exist
        os.makedirs(path, mode=0o755, exist_ok=False)
        os.chmod(path, 0o755)  # just in case mode above doesn't process
        print("The new directory is created! path:", path)

    # write file in sfs_output [1:] .csv
    arr = np.asarray(site_frequency_spectrum)
    pd.DataFrame(arr).to_csv("./sfs_output/sfs-output.csv")

    return site_frequency_spectrum


def frequency(sequence_dataframe, data):
    print("calculating frquencies")
    print(sequence_dataframe, "\n", data)

    for position in sequence_dataframe:
        print(sequence_dataframe[position].value_counts(normalize=True))

    return


def oberserved_heterozygosity():
    return


def inbreeding_coefficient(sequence_dataframe, data):
    print("calculating inbreeding coefficient")

    print(data)
    print(sequence_dataframe)
    individuals = data["IDV"].unique()
    heteroz_dict = {}

    for indv in individuals:
        col_to_concatenate = data["IDV"] == indv
        print(sequence_dataframe[col_to_concatenate])

        individual_dataframe = sequence_dataframe[col_to_concatenate]

        heteroz_dict[indv] = []

        for i, allele in enumerate(individual_dataframe):
            print(individual_dataframe[i][0], individual_dataframe[i][1])
            allele_combination = "".join(individual_dataframe[i].astype("str"))
            heteroz_dict[indv].append(allele_combination)

    diploid_df = pd.DataFrame.from_dict(heteroz_dict)
    print(diploid_df)

    for position in sequence_dataframe:
        frequencies = (
            pd.Series(sequence_dataframe[position])
            .value_counts(normalize=True)
            .tolist()
        )

        H_exp = 1
        for site_freq in frequencies:
            H_exp *= site_freq
        print(frequencies, H_exp)

    return


def nei_fst(sequence_dataframe, data, identifier):  # nei 1977
    print("caluculating nei_fst")
    # print(data)
    # print(sequence_dataframe)

    # identififer for subpopulaitons
    ID = identifier

    # initialize denom and numerator
    H_s_sub = {}
    H_t_pos = []

    # identify subpopulations
    subpopulations = data[ID].unique()

    # H_s
    for subpopulation in subpopulations:
        H_s_sub[subpopulation] = []  # positional array
        isSUBPOPULATION = data[ID] == subpopulation
        subpopulation_sequences = sequence_dataframe[isSUBPOPULATION]
        for i in range(sequence_dataframe.shape[1]):
            # print(i)
            # nomralize base freqs 1-(2pq+2pr+2qr+ 2ps+2qs+2rs
            base_position_frequency = subpopulation_sequences.iloc[:, i].value_counts(
                normalize=True
            )
            length = len(
                subpopulation_sequences.iloc[:, i].value_counts(normalize=True)
            )
            if length == 1:
                H_s_sub[subpopulation].append(0.001)
            else:
                H_s_sub[subpopulation].append(
                    expected_heterozygosity(base_position_frequency)
                )
    # print('Expected Heterozygosity for the each subpopulation', H_s_sub)    #for pos

    # H_t
    for i in range(sequence_dataframe.shape[1]):
        base_position_frequency = sequence_dataframe.iloc[:, i].value_counts(
            normalize=True
        )
        length = len(sequence_dataframe.iloc[:, i].value_counts(normalize=True))
        if length == 1:
            H_t_pos.append(0.001)
        else:
            H_t_pos.append(expected_heterozygosity(base_position_frequency))
    # print('Expected Heterozygosity for the whole population', H_t_pos)

    # F_st final
    F_st_pop = {}
    for subpopulation in subpopulations:
        F_st_pop[subpopulation] = []  # positional array
        for i, H_t in enumerate(H_t_pos):
            H_s = H_s_sub[subpopulation][i]

            # Ht-Hs / Ht
            numer = H_t - H_s
            denom = H_t
            F_st = numer / denom

            # 1 - Hs/Ht
            # F_st = 1
            # F_st -= H_s / H_t

            F_st_pop[subpopulation].append(F_st)

    N_k = len(subpopulations)  # number of demes
    sequence_range = range(len(F_st_pop[subpopulation]))

    # fig, ax = plt.subplots(N_k, figsize=(10, 3*N_k))
    # for i, subpopulation in enumerate(subpopulations):
    #    ax[i].scatter(y = F_st_pop[subpopulation], x = sequence_range)
    #     ax[i].set_xlabel("Position of deme" + str(subpopulation))
    #      ax[i].set_ylabel("F_st")
    #       ax[i].set_ylim(0,1)
    #
    #
    # plt.show()

    # fix fst avg over loci
    numerator = {}
    Fst_over_loci = {}
    for subpopulation in subpopulations:
        numerator[subpopulation] = []
        Fst_over_loci[subpopulation] = 0
        denomernator = []
        for i, H_t in enumerate(H_t_pos):
            H_s = H_s_sub[subpopulation][i]

            # Ht-Hs / Ht
            numer = H_t - H_s
            denom = H_t

            numerator[subpopulation].append(numer)
            denomernator.append(denom)
        num_avg = np.mean(numerator[subpopulation])
        den_avg = np.mean(denomernator)
        Fst_over_loci[subpopulation] = num_avg / den_avg

    for subpopulation, Fst in Fst_over_loci.items():
        print("Fst", subpopulation, ":", Fst)

    return


def expected_heterozygosity(normalized_value_counts):
    H_e = 1

    for i in range(len(normalized_value_counts)):
        H_e -= normalized_value_counts.iloc[i] * normalized_value_counts.iloc[i]

    return H_e


def Shannon_index(normalized_value_counts):
    print(normalized_value_counts)
    E_shannon = 1

    for i in range(len(normalized_value_counts)):
        E_shannon -= normalized_value_counts.iloc[i] * np.log(
            normalized_value_counts.iloc[i]
        )

    return E_shannon


def expected_heterozygosity_depricated(normalized_value_counts):
    print(normalized_value_counts)
    H_e = 0

    if len(normalized_value_counts) == 2:
        H_e += 2 * normalized_value_counts.iloc[0] * normalized_value_counts.iloc[1]
    else:
        for i in range(len(normalized_value_counts)):
            for j in range(len(normalized_value_counts)):
                if i < j:
                    H_e += (
                        2
                        * normalized_value_counts.iloc[i]
                        * normalized_value_counts.iloc[j]
                    )

    return H_e


def nei_chesser_pairwise_fst(sequence_dataframe, data, identifier):
    print("caluculating pairwise_fst, use (NC83 Eq. 15) and (NC83 Eq. 16)")
    # assumption demes in hardy-wieburg equi
    # nei and chesser
    # H_s = 2n_til / (2n_til - 1) * (1 - sigma_over_bases basefreqavg_overdemes*square)
    # H_s_i

    # harmonic mean on deme populations

    ##number of demes

    #


def wright_fis(sequence_dataframe, data, identifier):
    print("caluculating fis")
    # identififer for subpopulaitons
    ID = identifier

    # initialize denom and numerator
    H_s_sub = {}
    H_i_sub = {}

    # identify subpopulations
    subpopulations = data[ID].unique()

    # H_s
    for subpopulation in subpopulations:
        H_s_sub[subpopulation] = []  # positional array
        isSUBPOPULATION = data[ID] == subpopulation
        subpopulation_sequences = sequence_dataframe[isSUBPOPULATION]
        for i in range(sequence_dataframe.shape[1]):
            base_position_frequency = subpopulation_sequences.iloc[:, i].value_counts(
                normalize=True
            )
            length = len(
                subpopulation_sequences.iloc[:, i].value_counts(normalize=True)
            )
            if length == 1:
                H_s_sub[subpopulation].append(0.001)
            else:
                H_s_sub[subpopulation].append(
                    expected_heterozygosity(base_position_frequency)
                )

    # H_i
    for subpopulation in subpopulations:
        H_i_sub[subpopulation] = []  # positional array
        isSUBPOPULATION = data[ID] == subpopulation
        subpopulation_sequences = sequence_dataframe[isSUBPOPULATION]
        for i in range(sequence_dataframe.shape[1]):
            H_i_sub[subpopulation].append(0)
            position_data = pd.merge(
                subpopulation_sequences.iloc[:, i],
                data[isSUBPOPULATION],
                left_index=True,
                right_index=True,
            )
            position = position_data.columns.to_list()[0]
            grouped_df = position_data.reset_index().groupby(["IDV"])[position]

            for key, item in grouped_df:
                chrom_comb = grouped_df.get_group(key).index
                chr1_base = grouped_df.get_group(key)[chrom_comb[0]]
                chr2_base = grouped_df.get_group(key)[chrom_comb[1]]
                if chr1_base != chr2_base:
                    H_i_sub[subpopulation][i] += 1
            H_i_sub[subpopulation][i] /= sequence_dataframe.shape[1]
    # F_is final
    F_is_pop = {}
    for subpopulation in subpopulations:
        F_is_pop[subpopulation] = []  # positional array
        for i, H_i_deme in enumerate(H_i_sub):
            H_s = H_s_sub[subpopulation][i]
            H_i = H_i_sub[H_i_deme][i]
            # Ht-Hs / Ht
            numer = float(H_s) - float(H_i)
            denom = H_s
            F_is = numer / denom

            F_is_pop[subpopulation].append(F_is)

    N_k = len(subpopulations)  # number of demes
    sequence_range = range(len(F_is_pop[subpopulation]))

    # fig, ax = plt.subplots(N_k, figsize=(10, 3*N_k))
    # for i, subpopulation in enumerate(subpopulations):
    #    ax[i].scatter(y = F_is_pop[subpopulation], x = sequence_range)
    #    ax[i].set_xlabel("Position of deme" + str(subpopulation))
    #    ax[i].set_ylabel("F_is")
    #    ax[i].set_ylim(0,1)#
    #
    #
    #   plt.show()

    # fix fst avg over loci
    numerator = {}
    Fis_over_loci = {}
    denomernator = {}
    for subpopulation in subpopulations:
        numerator[subpopulation] = []
        Fis_over_loci[subpopulation] = 0
        denomernator[subpopulation] = []
        for i, H_i_deme in enumerate(H_i_sub):
            H_s = H_s_sub[subpopulation][i]
            H_i = H_i_sub[H_i_deme][i]
            # Ht-Hs / Ht
            numer = float(H_s) - float(H_i)
            denom = H_s
            F_is = numer / denom

            F_is_pop[subpopulation].append(F_is)
            numerator[subpopulation].append(numer)
            denomernator[subpopulation].append(denom)
        num_avg = np.mean(numerator[subpopulation])
        den_avg = np.mean(denomernator[subpopulation])
        Fis_over_loci[subpopulation] = num_avg / den_avg

    for i, subpopulation in enumerate(subpopulations):
        print("Fis", subpopulation, ":", np.mean(F_is_pop[subpopulation]))

    return


def by_deme_pairwise_fst(sequence_dataframe, data, identifier):
    print("caluculating pairwise_fst")

    # print(data)
    # print(sequence_dataframe)

    # identififer for subpopulaitons
    ID = identifier

    # initialize denom and numerator
    H_s_sub = {}
    H_ss_sub = {}
    H_t_pos = []

    # identify subpopulations
    subpopulations = data[ID].unique()

    # H_s
    for subpopulation in subpopulations:
        H_s_sub[subpopulation] = []  # positional array
        isSUBPOPULATION = data[ID] == subpopulation
        subpopulation_sequences = sequence_dataframe[isSUBPOPULATION]
        for i in range(sequence_dataframe.shape[1]):
            base_position_frequency = subpopulation_sequences.iloc[:, i].value_counts(
                normalize=True
            )
            length = len(
                subpopulation_sequences.iloc[:, i].value_counts(normalize=True)
            )
            if length == 1:
                H_s_sub[subpopulation].append(0.001)
            else:
                H_s_sub[subpopulation].append(
                    expected_heterozygosity(base_position_frequency)
                )

    # H_ss
    for subpopulation_i in subpopulations:
        for subpopulation_i_prime in subpopulations:
            comparison = subpopulation_i + "->" + subpopulation_i_prime
            if subpopulation_i != subpopulation_i_prime:
                H_ss_sub[comparison] = []
                isSUBPOPULATION_i = data[ID] == subpopulation_i
                isSUBPOPULATION_i_prime = data[ID] == subpopulation_i_prime
                isCOMPARISON = isSUBPOPULATION_i + isSUBPOPULATION_i_prime
                comparison_sequences = sequence_dataframe[isCOMPARISON]
                for i in range(sequence_dataframe.shape[1]):
                    base_position_frequencies = comparison_sequences.iloc[
                        :, i
                    ].value_counts(normalize=True)
                    H_ss_sub[comparison].append(
                        expected_heterozygosity(base_position_frequency)
                    )

    # H_t
    for i in range(sequence_dataframe.shape[1]):
        base_position_frequency = sequence_dataframe.iloc[:, i].value_counts(
            normalize=True
        )
        length = len(sequence_dataframe.iloc[:, i].value_counts(normalize=True))
        if length == 1:
            H_t_pos.append(0.001)
        else:
            H_t_pos.append(expected_heterozygosity(base_position_frequency))
    # print('Expected Heterozygosity for the whole population', H_t_pos)

    # F_st final
    F_st_pair = {}
    for subpopulation_i in subpopulations:
        for subpopulation_i_prime in subpopulations:
            if subpopulation_i != subpopulation_i_prime:
                comparison = subpopulation_i + "->" + subpopulation_i_prime
                F_st_pair[comparison] = []  # positional array
                for i, H_t in enumerate(H_ss_sub[comparison]):
                    H_s = H_s_sub[subpopulation_i][i]
                    F_st = 1
                    F_st -= H_s / H_t
                    F_st_pair[comparison].append(F_st)

    # print("F_st array", F_st_pair)
    # print('add Fst plots per deme/subpop/subspec')
    # length unique demes N_k
    # subplots(N_k,figsize=(10,3*N_k))
    N_k = len(subpopulations)  # number of demes
    sequence_range = range(len(F_st_pair[comparison]))
    num_axes = math.comb(N_k, 2) * 2
    # fig, ax = plt.subplots(num_axes, figsize=(10, 3 * num_axes))
    # onAX = 0
    # for  i, subpopulation_i in enumerate(subpopulations):
    #    for j, subpopulation_i_prime in enumerate(subpopulations):
    #        comparison = subpopulation_i + '->' + subpopulation_i_prime
    #        if subpopulation_i != subpopulation_i_prime:
    #            ax[onAX].scatter(y = F_st_pair[comparison], x = sequence_range)
    #            ax[onAX].set_xlabel("Position of deme " + str(comparison))
    #            ax[onAX].set_ylabel("F_st")
    #            ax[onAX].set_ylim(0,1)
    #            onAX += 1
    # plt.show()

    # fix fst avg over loci
    numerator = {}
    Fst_over_loci = {}
    denomernator = {}
    for subpopulation_i in subpopulations:
        for subpopulation_i_prime in subpopulations:
            if subpopulation_i != subpopulation_i_prime:
                comparison = subpopulation_i + "->" + subpopulation_i_prime
                numerator[comparison] = []
                Fst_over_loci[comparison] = 0
                denomernator[comparison] = []
                for i, H_t in enumerate(H_t_pos):
                    H_s = H_s_sub[subpopulation][i]

                    # Ht-Hs / Ht
                    numer = H_t - H_s
                    denom = H_t

                    numerator[comparison].append(numer)
                    denomernator[comparison].append(denom)
                num_avg = np.mean(numerator[comparison])
                den_avg = np.mean(denomernator[comparison])
                Fst_over_loci[comparison] = num_avg / den_avg

    for comparison, Fst in Fst_over_loci.items():
        print("Fst", comparison, ":", Fst)

    return


def weir_goudet_population_specific_fst(sequence_dataframe, data, identifier):
    print("caluculating weir_goudet_population_specific_fst")
    # print(data)
    # print(sequence_dataframe)

    # identififer for subpopulaitons
    ID = identifier

    # initialize denom and numerator
    M_w_sub = {}
    M_B_pos = []

    # identify subpopulations
    subpopulations = data[ID].unique()

    # M_W
    for subpopulation in subpopulations:
        M_w_sub[subpopulation] = []  # positional array of seq
        isSUBPOPULATION = data[ID] == subpopulation
        deme_population_size = isSUBPOPULATION.sum()  # eta
        subpopulation_sequences = sequence_dataframe[isSUBPOPULATION]
        M_w_constanta = (2 * deme_population_size) / (2 * deme_population_size - 1)
        M_w_constantb = (1) / (2 * deme_population_size - 1)
        for i in range(sequence_dataframe.shape[1]):
            base_position_frequency = subpopulation_sequences.iloc[:, i].value_counts(
                normalize=True
            )
            square_base_frequency = np.square(base_position_frequency)
            sum_square_base_frequency = square_base_frequency.sum()
            M_w_sub[subpopulation].append(
                M_w_constanta * (sum_square_base_frequency - M_w_constantb)
            )
    # print('M_W',M_w_sub)
    # M_B
    M_B_constant = 1 / (len(subpopulations) * (len(subpopulations) - 1))

    for i in range(sequence_dataframe.shape[1]):
        M_B = 0

        for subpopulation_i in subpopulations:
            for subpopulation_i_prime in subpopulations:
                if subpopulation_i != subpopulation_i_prime:
                    isSUBPOPULATION_i = data[ID] == subpopulation_i
                    isSUBPOPULATION_i_prime = data[ID] == subpopulation_i_prime
                    subpopulation_i_sequences = sequence_dataframe[isSUBPOPULATION_i]
                    subpopulation_i_prime_sequences = sequence_dataframe[
                        isSUBPOPULATION_i_prime
                    ]

                    base_position_frequency_i = subpopulation_i_sequences.iloc[
                        :, i
                    ].value_counts(normalize=True)
                    base_position_frequency_i_prime = (
                        subpopulation_i_prime_sequences.iloc[:, i].value_counts(
                            normalize=True
                        )
                    )
                    base_position_i_prime = subpopulation_i_prime_sequences.iloc[:, i]

                    # sum products of base frequency
                    bases_in_i = base_position_frequency_i.index
                    bases_in_i_prime = base_position_frequency_i_prime.index
                    for j in np.intersect1d(bases_in_i, bases_in_i_prime):
                        M_B += (
                            base_position_frequency_i[j]
                            * base_position_frequency_i_prime[j]
                        )
        M_B_pos.append(M_B * M_B_constant)

    for deme, M_w_i_pos in M_w_sub.items():
        numerator = 0
        denomenator = 0
        seq_length = len(M_w_i_pos)
        for l in range(seq_length):
            numerator += M_w_i_pos[l] - M_B_pos[l]
            denomenator += 1 - M_B_pos[l]

        numerator /= seq_length
        denomenator /= seq_length
        print("psFst", deme, ":", numerator / denomenator)
    return numerator / denomenator


def decrepit_code():
    for subpopulation_i in subpopulations:
        for subpopulation_i_prime in subpopulations:
            if subpopulation_i != subpopulation_i_prime:
                isSUBPOPULATION_i = data[ID] == subpopulation_i
                isSUBPOPULATION_i_prime = data[ID] == subpopulation_i_prime
                subpopulation_i_sequences = sequence_dataframe[isSUBPOPULATION_i]
                subpopulation_i_prime_sequences = sequence_dataframe[
                    isSUBPOPULATION_i_prime
                ]
                for i in range(sequence_dataframe.shape[1]):
                    M_B = 0
                    base_position_frequency_i = subpopulation_i_sequences.iloc[
                        :, i
                    ].value_counts(normalize=True)
                    base_position_frequency_i_prime = (
                        subpopulation_i_prime_sequences.iloc[:, i].value_counts(
                            normalize=True
                        )
                    )
                    base_position_i_prime = subpopulation_i_prime_sequences.iloc[:, i]

                    # sum products of base frequency
                    bases_in_i = base_position_frequency_i.index
                    bases_in_i_prime = base_position_frequency_i_prime.index
                    for j in np.intersect1d(bases_in_i, bases_in_i_prime):
                        M_B += (
                            base_position_frequency_i[j]
                            * base_position_frequency_i_prime[j]
                        )
                    M_B_pos.append(M_B * M_B_constant)
    print("M_B", M_B_pos)
