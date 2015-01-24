#!/user/bin/env python
"""Filters a pVAAST report by a set frequency cutoff using EXaC data, currently
 only works with SNVs."""

import csv, argparse

__author__="Brett J. Kennedy"
__date__="January 24 2015"

def parse_args():
    """parses the args"""
    pass

def exac_freq(chrom,pos,allele,exac):
    """finds the frequency of an allele in the EXaD database"""
    pass

def genotype_freq(info,exac):
    """combines the exac frequency of alleles in the info field to 
    calcualte a genotype frequency for the scored allels"""
    pass

def write_out(out,line):
    """writes the to the output"""
    pass

def rerank(genes):
    """reranks the genes after those which fail the EXaC filter
    are removed"""
    pass

def main()