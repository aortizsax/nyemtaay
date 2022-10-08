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
from nyemtaay.calculate import populationgeneticstats
from nyemtaay.mathlib import sterling
from nyemtaay.tests.nuetrality import tajimas_d

def main():
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument(
        "-f",
        "--fastafiles",
        action="store",
        nargs="+",
        metavar="FILE",
        help="Path to source file(s).",
    )
    parser.add_argument(
        "-m",
        "--metadata",
        action="store",
        nargs="+",
        metavar="FILE",
        help="Path to source file(s).",
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
        help="Prefix for output files [default=%(default)s].",
    )
    parser.add_argument(
        "-calc",
        "--calculate",
        action="store",
        default="False",
        help="Prefix for output files [default=%(default)s].",
    )

    args = parser.parse_args()
    print("Parsing")
    
    # use parser modules
    # pass list of fasta files to fasta parser
    sequence_matrix = read_fasta_files(args.fastafiles)
    
    # pass metadata to its parser
    data_matrix = read_metadata(args.metadata[0],args.header)
    print("Done parsing")
    
    # convert matrix to dataframe with indexes matching metadata
    sequence_dataframe = to_dataframe(sequence_matrix, data_matrix)
    
    sequence_dataframe.wright_fst = populationgeneticstats.wright_fst(sequence_dataframe, data_matrix)
    
    sequence_dataframe.segregating_sites  = populationgeneticstats.number_segregating_sites(sequence_dataframe,data_matrix)
    
    
    tajimas_d(sequence_dataframe,data_matrix)
    
    #sequence_dataframe.nuc_div = populationgeneticstats.nucleotide_diversity(sequence_dataframe,data_matrix)

        
    
    sequence_dataframe.sfs = populationgeneticstats.frequency_spectrum(sequence_dataframe,data_matrix)


    #populationgeneticstats.inbreeding_coefficient(sequence_matrix,data_matrix)


if __name__ == "__main__":
    main()
