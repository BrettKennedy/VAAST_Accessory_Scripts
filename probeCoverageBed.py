import csv
import argparse
import os
import numpy as np
import pandas as pd
import operator


def get_args():
    parser = argparse.ArgumentParser(description="takes a list of gatk "
                                     "\"sample_interval_summary\" files and "
                                     "generates a bed file sorted by the mean"
                                     " depth of coverage across all samples, "
                                     "also contains columns for the coverage "
                                     "at that position for each sample")
    parser.add_argument("summaries", nargs='+', help="space seperated list"
                        " of sample_interval_summary files")
    parser.add_argument("--output", help="name of output bed file")
    return parser.parse_args()


def parseSummary(summary, coverageDict):
    """Parses the GATK generated sample_interval_summary into """
    if coverageDict.keys():
        with open(summary) as t:
            for line in csv.reader(t, delimiter="\t"):
                if line[0] == "Target":
                    continue
                try:
                    coverageDict[line[0]].append(float(line[2]))
                except KeyError:
                    print "Coverage files not created with the same bed"
    else:
        with open(summary) as t:
            for line in csv.reader(t, delimiter="\t"):
                if line[0] == "Target":
                    continue
                coverageDict[line[0]] = [float(line[2])]
    return 0


def generateHeader(summaries, output):
    """ Writes out the header file by column to the output bed file """
    header = ["#chr", "region_start", "region_end", "mean_coverage"]
    header.extend([i.strip("_newBed.sample_interval_summary").strip("\/")
                   for i in summaries])
    header.append("\n")
    output.write("\t".join(header))
    return 0


def generateBed(summaries, output):
    """ Generates the bed files """
    coverageDict = {}
    meanDict = {}
    for summary in summaries:
        parseSummary(summary, coverageDict)
    for region in coverageDict.keys():
        meanDict[region] = np.mean(coverageDict[region])
    with open(output, 'w') as out:
        generateHeader(summaries, out)
        for region in sorted(meanDict, key=lambda k: meanDict[k]):
            covReg = [str(i) for i in coverageDict[region]]
            reg = "%s\t%s" % (region.split(':')[0],
                              "\t".join(region.split(':')[1].split('-')))
            o = "%s\t%1.3f\t%s\n" % (reg, meanDict[region], "\t".join(covReg))
            out.write(o)


def main():
    args = get_args()
    generateBed(args.summaries, args.output)


if __name__ == "__main__":
    main()
