#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 09:16:34 2019

@author: meike
"""

#--outdir /home/meiker/git/prokka_annotation/streptococcus_02451
# 'streptococcus_11895', 'streptococcus_11939']
# 'lactococcus_00188'
#wget ftp://ftp.patricbrc.org/genomes/1739284.3/1739284.3.fna -O /home/meiker/git/genomes/streptocuccus_11846.fna

 x_min = min(x_in)
    x_max = max(x_in)
    x = [[] for _ in range(n_bins)]
    for a in x_in:
        # compute the bin number for value a
        n = int(float(a - x_min) / (x_max - x_min + 1.0) * n_bins)
        x[n].append(a)
    return x  # x is a binned list of elements from x_in