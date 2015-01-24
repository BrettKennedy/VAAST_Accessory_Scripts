#!/user/bin/env python
"""Filters a pVAAST report by a set frequency cutoff using EXaC data, currently
 only works with SNVs."""

import csv, argparse, operator
import tabix

__author__="Brett J. Kennedy"
__date__="January 24 2015"

def parse_args():
  """parses the args"""
  parser=argparse.ArgumentParser(description="Filters a pVAAST report by a set\
   frequency cutoff using EXaC data")
  parser.add_argument("pvaast",help="pVAAST report to be filtered")
  parser.add_argument("output",help="output for filtered and reranked report")
  parser.add_argument("cutoff",help="genotype frequency cutoff",type=float)
  parser.add_argument("exac",help="location tabix indexed EXaC database")
  return parser.parse_args()

def exac_freq(chrom,pos,allele,exac):
  """finds the frequency of an allele in the EXaD database"""
  freq=0
  try:  equer=exac.query(chrom,pos-1,pos)
  except: return freq
  for line in equer:
    ai=0
    alleles=line[4].split(',')
    if len(alleles)>=1:
      for a in range(len(alleles)):
        if alleles[a]==allele:
          ai=a
    if int(line[1])!=pos:
      continue
    info=line[7].split(';')
    count=float(info[0].split('=')[1].split(',')[ai])
    total=float(info[12].split('=')[1])
  freq=count/total
  return freq

def genotype_freq(info,exac):
  """combines the exac frequency of alleles in the info field to 
  calcualte a genotype frequency for the scored allels"""
  freqs=[]
  for var in info:
    posvar=var.split(';')[0]
    allele=var.split(';')[-3].split('->')[1]
    chrom,pos=posvar.strip('chr').split(':')
    freqs.append(exac_freq(chrom,pos,allele,exac))
  if len(freqs)==1:
    return freqs[0]**2
  elif len(freqs)==2:
    return freqs[0]*freqs[1]
  elif len(freqs)==3:
    return freqs[0]*freqs[1]*freqs[2]
  elif len(freqs)>=4:
    print "4+ alleles? thats unpossible"
    return 0

def write_out(output,ranked_list):
  """writes the to the output"""
  out=open(output,'w')
  out.write("RANK\tGene\tp-value\tp-value-ci\tScore\tLOD\tVariants\n")
  rank=0
  for line in ranked_list:
    rank+=1
    o=[rank]+line
    o=[str(i) for i in o]
    o="\t".join(o)+"\n"
    out.write(o)
  out.close()

def rerank(genes):
  """reranks the genes after those which fail the EXaC filter
  are removed"""
  list_to_sort=[]
  for g in genes:
    list_to_sort.append(genes[g])
  lis=sorted(list_to_sort, key=operator.itemgetter(1))
  return lis

def main(args):
  genes={}
  exac=tabix.Tabix(args.exac)
  with open(args.pvaast) as t:
    for line in csv.reader(t,delimiter="\t"):
      if line[0]=="RANK": continue
      line[2]=float(line[2])
      if line[2]=="1":
        genes[line[1]]=line[1:]
        continue
      posinfo=line[6:]
      gt_freq=genotype_freq(posinfo,exac)
      print gt_freq
      if gt_freq<args.cutoff:
        genes[line[1]]=line[1:]
      else: continue
  reranked_list=rerank(genes)
  write_out(args.output,reranked_list)

if __name__=="__main__":
  args=parse_args()
  main(args)