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


def end_pt_finder(genetic_map_filename, number_of_chromosomes=10):
    genetic_map = pd.read_csv(genetic_map_filename, sep='\t', index_col=None)
    end_pts = [0]
    for i in range(1, len(genetic_map)):
        if genetic_map[i-1][0] != genetic_map[i][0]:
            end_pts.append(i)
        else:
            pass
    end_pts.append(len(genetic_map))
    chrom_lengths = [(end_pts[i] - end_pts[i-1]) for i in range(1, number_of_chromosomes+1)]
    del end_pts[0]
    return end_pts, chrom_lengths

# Both qtl and ae file readers/writers need serious reconstruction.


def ae_file_reader(pop, filename):
    ae_dict = col.OrderedDict()
    with open(filename, newline='\n') as ae_file:
        ae_reader = csv.reader(ae_file, delimiter=',', quotechar='"')
        for row in ae_reader:
            splitrows = row[0].split(",")
            stripped_sublocus = splitrows[0].lstrip('(')
            stripped_nucleotide = splitrows[1].rstrip(')')
            sl = int(stripped_sublocus)
            nt = int(stripped_nucleotide)
            effect = float(row[1])
            ae_dict[sl, nt] = effect
    pop.dvars().alleleEffects = ae_dict


def parse_recombination_rates(genetic_map_filename):
    """
    Returns a list of crossover probabilities from a genetic map measured in centimorgans.
    """
    genetic_map = pd.read_csv(genetic_map_filename, sep='\t', index_col=None)
    genetic_map.drop(['locus', 'agpv2', 'namZmPRDA', 'namZmPRDS'], axis=1, inplace=True)
    genetic_map = np.array(genetic_map)
    recombination_rates = col.OrderedDict()
    for i in range(1, len(genetic_map), 1):
        if genetic_map[i-1][0] == genetic_map[i][0]:
            recombination_rates[i] = np.divide(np.abs(genetic_map[i][1] - genetic_map[i-1][1]), 100)
        elif genetic_map[i-1][0] != genetic_map[i][0]:
            recombination_rates[i] = 0.0
    recombination_rates[len(genetic_map)] = 0.0
    return list(recombination_rates.values())


def parse_genotype_matrix(genotype_matrix_filename: str, columns_to_drop='popdata', index_of_first_nonfounder=105):
    genotype_matrix = pd.read_csv(genotype_matrix_filename, sep='\t', index_col=0, low_memory=False)
    droppable_individuals = list(genotype_matrix.index[index_of_first_nonfounder:])
    genotype_matrix = genotype_matrix.drop(droppable_individuals, axis=0)
    if columns_to_drop is not None:
        genotype_matrix = genotype_matrix.drop(columns_to_drop, axis=1)
    return genotype_matrix
