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
# library for tests of neutrality
#
###########################################################################
## Reader Implementation

# Libraries
import numpy as np
import matplotlib.pyplot as plt
from nyemtaay.calculate.populationgeneticstats import (
    nucleotide_diversity,
    number_segregating_sites,
)
from nyemtaay.mathlib import sterling


def tajimas_d(sequence_dataframe, data):
    # Tajimas 1989
    # test with beta distribution
    print("Testing tajimas_D")
    D_numerator = nucleotide_diversity(sequence_dataframe, data)
    number_seg_sites = number_segregating_sites(sequence_dataframe, data)
    n = len(data.index)
    sterling_number = sterling.first_order_inverse(n)

    D_numerator -= number_seg_sites / sterling_number

    D_denominator = 1
    # e_1(S) + e_2(S)(S-1)
    # e_1 = 1/a_1((n+1)/(3(n-1))-(1/a_1))
    e_1 = n + 1
    e_1 /= 3 * (n - 1)
    e_1 -= 1 / sterling.first_order_inverse(n)
    e_1 *= 1 / sterling.first_order_inverse(n)

    # e_2 = (1/(a_1+a_2)) ((2(n^2+n+3)/(9n(n-1))) - ((n+2)/(na_1)) + (a_1/a_2))
    # e_2 = a (b -c +d)
    # a = (1/(a_1+a_2))
    a = 1 / (sterling.first_order_inverse(n) + sterling.second_order_inverse(n))
    # b = (2(n^2+n+3)/(9n(n-1)))
    b = 2 * (n ^ 2 + n + 3) / (9 * n * (n - 1))
    # c = ((n+2)/(na_1))
    c = (n + 2) / (n * sterling.first_order_inverse(n))
    # d = (a_1/a_2)
    d = sterling.first_order_inverse(n) / sterling.second_order_inverse(n)

    e_2 = a * (b - c + d)

    D_denominator = e_1 * number_seg_sites + e_2 * number_seg_sites * (
        number_seg_sites - 1
    )

    tajimas = D_numerator / D_denominator

    print("Tajimas numerator", tajimas)

    return tajimas
