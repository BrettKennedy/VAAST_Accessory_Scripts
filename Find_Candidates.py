#!/usr/bin/env python
"""Finds candidate genes in a pVAAST (or a phevor2) report above a given rank cutoff\
returns variant info"""

__author__="Brett J. Kennedy"
__date__="January 30, 2015"

import csv, argparse, os, sys
from exac_sieve import exac_freq
from exac_sieve import genotype_freq
import tabix

def parse_args():
  """parses the arguments"""
  parser=argparse.ArgumentParser(description="Finds candidate genes in a \
phevor or pVAAST report above a given rank cutoff returns variant info")
  parser.add_argument("pvaast",help="pVAAST simple report")
  parser.add_argument("cutoff",help="rank cutoff below which genes cannot be reported",type=int)
  parser.add_argument("candidates",help="flat text file of candidate genes")
  parser.add_argument("--phevor",help="if included rank cutoff will be based on phevor rank",default=None)
  parser.add_argument("--exac",help="location of exac DB, if included exac GT frequencies are reported",default=None)
  return parser.parse_args()

def write_header(args):
  header=["Sample ID",\
  "Chromosome position",\
  "Function", "Gene",\
  "Protien change",\
  "Protein function",\
  "Polyphen",\
  "VAAST Score",\
  "VAAST p-value",\
  "pVAAST rank",]
  if args.exac:
    header.append("ExAC GT Frequency")
  if args.phevor:
    header.append("Phevor2 rank")
  o="\t".join(header)+"\n"
  sys.stdout.write(o)

def outputter(genesD,args):
  for gene in genesD:
    o=[str(i) for i in genesD[gene]]
    o="\t".join(o)+"\n"
    sys.stdout.write(o)

def parse_cands(cand_list):
  list=[]
  with open(cand_list) as t:
    for line in csv.reader(t):
      list.append(line[0])
  return list

def parse_pvaast(args,cands,exac=None):
  pvaastD={}
  cutoff=args.cutoff
  if args.phevor:
    cutoff=26000
  with open(args.pvaast) as t:
    for line in csv.reader(t,delimiter="\t"):
      if line[0]=="RANK":  continue
      if line[1] not in cands:	continue
      if int(line[0])>=cutoff:  break
      if float(line[2])==1:  break
      if float(line[5])==0:  continue
      id=args.pvaast.split("_recessive_")[0]	
      posinfo=" ".join(get_pos_info(line[6:]))
      gene=line[1]
      if posinfo==[]:  continue
      pvaastD[gene]=[id,posinfo,"",gene,"","","",line[4],line[2],line[0]]
      if args.exac:
        gt_freq=genotype_freq(line[6:],exac)
        pvaastD[gene].append(gt_freq)
  return pvaastD
      
def get_pos_info(info):
  posinfo=[]
  for position in info:
    pinfo=position.split(';')[0].strip('chr').split(':')
    ntinfo=position.split(';')[3].split('->')
    pinfo.extend(ntinfo)
    pinfo=" ".join(pinfo)
    posinfo.append(pinfo)
  return posinfo

def parse_phevor(args,pvgenes):
  phdict={}
  with open(args.phevor) as t:
    for line in csv.reader(t,delimiter="\t"):
      if "#" in line[0]:  continue
      line=[i.strip() for i in line]
      rank=int(int(line[0])+1)
      if rank > args.cutoff:  break 
      gene=line[1]
      if gene in pvgenes.keys():
        phline=pvgenes[gene]
        phline.append(rank)
        phdict[gene]=phline
  if phdict!={}:
    outputter(phdict,args)

def main():
  args=parse_args()
  cand_list=parse_cands(args.candidates)
  #write_header(args)
  if args.exac:
    exac=tabix.Tabix(args.exac)
    pvgenes=parse_pvaast(args,cand_list,exac)
  else:
    pvgenes=parse_pvaast(args,cand_list)
  if args.phevor:
    parse_phevor(args,pvgenes)
  else:
    outputter(pvgenes,args)

if __name__=="__main__":
  main()
