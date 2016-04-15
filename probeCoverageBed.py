import csv
import argparse
import os
import numpy as np
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser(description="takes a list of gatk "
                                     "\"sample_interval_summary\" files and "
                                     "generates a bed file sorted by the mean "
                                     "depth of coverage across all samples, "
                                     "also contains columns for the coverage "
                                     "at that position for each sample")
    parser.add_argument("summaries", nargs='+', help="comma seperated list"
                        " of sample_interval_summary files")
    parser.add_argument("output", help="name of output bed file")
    return parser.parse_args()


def parseSummary(summary, coverageFrame):
    """Parses the GATK generated sample_interval_summary into """
    if coverageFrame:
        with open(summary) as t:
            for line in csv.reader(t, delimiter="\t"):

    else:
        with open(summary) as t:
            for line in csv.reader(t, delimiter="\t")

    pass


def generateHeader(summaries, output):
    """ Writes out the header file by column to the output bed file """
    pass


def generateBed(summaries, output):
    """ Generates the bed files """
    pass


def main():
    args = get_args()
    generateBed(args.summaries, args.output)


if __name__ == "__main__":
    main()
