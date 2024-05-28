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

import os
import pathlib
import sys
import argparse

# import modules front end parser
from nyemtaay.parse.parser import read_fasta_files, read_metadata, to_dataframe
from nyemtaay.calculate import populationgeneticstats, informationtheory
from nyemtaay.mathlib import sterling
from nyemtaay.tests.nuetrality import tajimas_d


def main():
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument(
        "-f",
        "--fastafiles",
        action="store",
        nargs="+",
        required=True,
        metavar="FILE",
        help="Path to source file(s).",
    )
    parser.add_argument(
        "-m",
        "--metadata",
        action="store",
        nargs="+",
        required=True,
        metavar="FILE",
        help="Path to metadata file(s).",
    )
    parser.add_argument(
        "-s",
        "--subpopulation-identifier",
        action="store",
        required=True,
        default="String",
        help="Column in the metadata that you wish to use as deme ID for analysis (sampling site deme, sampling tissue, etc)",
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        action="store",
        default="output",
        help="Prefix for output files [default=%(default)s].",
    )
    parser.add_argument(
        "-hd",
        "--header",
        action="store",
        default=0,
        help="int for header line in metadata [default=%(default)s].",
    )
    parser.add_argument(
        "-calc",
        "--calculate",
        action="store",
        nargs="+",
        metavar="FUNCTION",
        help="List of analysis to perform. View readme on github.com/aortizsax/nyemtaay #move the following to readme'Shannon_index', 'pairwise_fst','expected_heterozygosity', 'frequency_spectrum', 'inbreeding_coefficient','nei_chesser_pairwise_fst', 'nei_fst','nucleotide_diversity', 'number_segregating_sites', 'oberserved_heterozygosity', 'segregating_sites', 'weir_goudet_population_specific_fst', 'wright_fis'. A by deme or a by deme and allele can be performed on 'shannon_entropy', 'directionality analysis', 'jsd', normalized_jsd, [default=%(default)s].",
    )

    parser.add_argument(
        "-dm",
        "--display-map",
        action="store_true",
        default=False,
        help="Map display option. With lat and long rows in the metadata, Nyemtaay will overlay deme anaylsis results to a map.",
    )

    #print(dir(populationgeneticstats), dir(informationtheory))

    args = parser.parse_args()
    print("Parsing")
    #print(dir(args))


    # Parse
    # use parser modules
    # pass list of fasta files to fasta parser
    sequence_matrix = read_fasta_files(args.fastafiles)

    # pass metadata to its parser
    data_matrix = read_metadata(args.metadata[0], args.header)
    print("Done parsing")

    # convert matrix to dataframe with indexes matching metadata
    sequence_dataframe = to_dataframe(sequence_matrix, data_matrix)
    print(sequence_dataframe)


##    #Analysis
##    num_segregating_sites = populationgeneticstats.number_segregating_sites(
##                            sequence_dataframe, data_matrix)
##
##    tajimas_D = tajimas_d(sequence_dataframe, data_matrix)
##
##    pi_nuc_div = populationgeneticstats.nucleotide_diversity(
##        sequence_dataframe, data_matrix
##    )
##
###    sfs = populationgeneticstats.frequency_spectrum(
###        sequence_dataframe, data_matrix
###    )
##
##    populationgeneticstats.nei_fst(
##        sequence_dataframe, data_matrix, args.subpopulation_identifier
##    )
    # populationgeneticstats.inbreeding_coefficient(sequence_matrix,data_matrix)

    #####################################################################

    #add map display
    populationgeneticstats.nei_fst(sequence_dataframe, data_matrix, args.subpopulation_identifier)
    populationgeneticstats.by_deme_pairwise_fst(sequence_dataframe, data_matrix, args.subpopulation_identifier)

    (gamete_probabilities, population_dict,pos_allele_probablity_deme_dict) = informationtheory.sequences_to_gamete_prob(sequence_dataframe, data_matrix, args.subpopulation_identifier, 'complete graph')
    geo = 'complete graph'

##    #by gamete
    (shannon_entropy_by_deme, weighted_bool) = informationtheory.demes_shannon_entropy(gamete_probabilities,geo)


    #add map display
    (norm_jsd_by_comp, perc_network) = informationtheory.demes_norm_jsd(gamete_probabilities, data_matrix, population_dict,geo)

    #add map display
    (I_by_comp , w_dict) = informationtheory.demes_info_flow_direction(gamete_probabilities, population_dict,perc_network,geo)
###
    (gamete_probabilities, population_dict) = informationtheory.sequences_to_random_deme_combinations(sequence_dataframe, data_matrix, args.subpopulation_identifier,perc_network,'complete graph')
    #add map display
    (IR_dict_by_comparison, weighted_bool) = informationtheory.randomized_information_flow_directionality(gamete_probabilities, I_by_comp,geo)

    print(I_by_comp, IR_dict_by_comparison, w_dict)


#    networkx_edge_list = informationtheory.information_theory_dict_to_pd_df(IR_dict_by_comparison,weighted_bool,geo)
#    informationtheory.plot_bidirectional_metric(networkx_edge_list)

    #by posisitional allele
    #informationtheory.demes_allele_shannon_entropy(pos_allele_probablity_deme_dict)
#    informationtheory.demes_allele_jsd(pos_allele_probablity_deme_dict, population_dict)
#    informationtheory.demes_allele_norm_jsd(pos_allele_probablity_deme_dict, population_dict)
#   (I_dict , w_dict) informationtheory.demes_allele_info_flow_direction(pos_allele_probablity_deme_dict, population_dict)



if __name__ == "__main__":
    main()
