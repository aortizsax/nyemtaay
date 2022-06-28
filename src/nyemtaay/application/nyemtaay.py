#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
## Copyright (c) 2022 Adrian Ortiz.
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
from nyemtaay.parse.parse_io import read_fasta_files, read_population_data,polymorphism_read_csv


def main():
    parser = argparse.ArgumentParser(description=None)
    # parser.add_argument(
        # "-sp",
        # "--sequence_polymorphism_files",
        # action="store",
        # nargs="+",
        # metavar="FILE",
        # help="Path to source file(s).",
    # )
    parser.add_argument(
        "-p",
        "--population_data",
        action="store",
        nargs="+",
        metavar="FILE",
        help="Path to population structure file(s).",
    )
    parser.add_argument(
        "-f",
        "--fasta_files",
        nargs='+',
        default="None",
        help="Prefix for output files [default=%(default)s].",
    )
    parser.add_argument(
        "-sp",
        "--sequence_polymorphism_files",
        action="store",
        default="None",
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
        default="0",
        help="header line for seq polym files [default=%(default)s].",
    )

    args = parser.parse_args()
    print("Parsing")
    if not isinstance(args.fasta_files,type(None)): 
        print('fasta_files')
        print(args.fasta_files,args.header)
        # use parser modules
        # pass list of fasta files to fasta parser
        read_fasta_files(args.fasta_files)
    

    else: 
        print('sp')
        print(args.sequence_polymorphism_files,args.header)
        # use parser modules
        # pass list of fasta files to fasta parser
        
        polymorphism_read_csv(args.sequence_polymorphism_files[0],int(args.header))
    

    if args.population_data:
        print('population')
        # pass population_data to its parser
        read_population_data(args.population_data)
    print("Done parsing")


if __name__ == "__main__":
    main()
