#! usr/bin/python
__author__ = 'John J. Dougherty III'
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import csv
import collections as col


def genotype_parser(genotype_matrix_filename, columns_to_drop=None,
                    droppable_individuals=None):
    """
    Specific to raw NAM files.
    Unsure if this is worth keeping here.
    However, I am unsure of how to implement a flexible parser module without assuming that the file has a specific
    structure already.
    Obtains marker data from a file which tabular data file the 26 inbred lines
    Assumes that genotype matrix will be in a format of Individual x Genotype Format
        Locus1  Locus2  Locus3  ...
    Ind1
    Ind2
    Ind3
    ...
    If there are any droppable columns/individuals specify that in the args.
    Low memory set to
    """
    genotype_matrix = pd.read_csv(genotype_matrix_filename, sep='\t', index_col=0, low_memory=False)
    if droppable_individuals is not None:
        genotype_matrix = genotype_matrix.drop(droppable_individuals, axis=0)
    if columns_to_drop is not None:
        genotype_matrix = genotype_matrix.drop(columns_to_drop, axis=1)
    return genotype_matrix

def parse_genotype_matrix(genotype_matrix_filename: str, columns_to_drop='popdata', index_of_first_nonfounder=105):
    genotype_matrix = pd.read_csv(genotype_matrix_filename, sep='\t', index_col=0, low_memory=False)
    droppable_individuals = list(genotype_matrix.index[index_of_first_nonfounder:])
    genotype_matrix = genotype_matrix.drop(droppable_individuals, axis=0)
    if columns_to_drop is not None:
        genotype_matrix = genotype_matrix.drop(columns_to_drop, axis=1)
    return genotype_matrix
