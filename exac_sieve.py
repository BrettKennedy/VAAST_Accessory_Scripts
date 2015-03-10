#!/user/bin/env python
"""Filters a pVAAST report by a set frequency cutoff using EXaC data, currently
 only works with SNVs."""

import csv, argparse, operator
import tabix

__author__="Brett J. Kennedy"
__date__="January 24 2015"

def parse_args():
  """parses the args"""
  parser=argparse.ArgumentParser(description="Filters a trio pVAAST report by a\
   set frequency cutoff using EXaC data.  Currently, does not work for pVAAST runs\
   that include indels.  Uses expected genotype frequencies for recessives.")
  parser.add_argument("pvaast",help="pVAAST report to be filtered")
  parser.add_argument("output",help="output for filtered and reranked report")
  parser.add_argument("cutoff",help="genotype frequency cutoff",type=float)
  parser.add_argument("exac",help="location tabix indexed EXaC database")
  parser.add_argument("--Filter_Hets",help="filters simple heterozygous genotypes (not transhets)",dest="hets",action='store_true',default=False)
  pheno=parser.add_mutually_exclusive_group()
  pheno.add_argument("--phevor",help="optionally, pass a phevor report along\
   with the pVAAST report, output will now be a phevor report",default=None)
  pheno.add_argument("--grouper",help="pass a grouper report with the pVAAST\
   report, output will now be ExAC filtered grouper output",default=None)
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
    for field in info:
      try:
        idx=field.split('=')[0]	
        value=field.split('=')[1]
      except:  continue
      if idx=='AC_Adj':
        count=float(value.split(',')[ai])
      if idx=='AN_Adj':
        total=float(value)
    if total==0:
      return 0
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
    freq=exac_freq(chrom,int(pos),allele,exac)
    if freq>0:
      freqs.append(freq)
  if freqs==[]:
    return 0.0
  elif chrom=="X":
    return freqs[0]/2
  elif len(freqs)==1:
    tarcount=info[0].split(',')[-1]
    if tarcount=="1":
      return freqs[0]
    if tarcount=="2":
      return freqs[0]**2
  elif len(freqs)>=2:
    return freqs[0]*freqs[1]

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

def grouper_in_out(grouper,genes,output):
  out=open(output,'w')
  ranker=0
  with open(grouper) as t:
    out.write(t.readline())
    out.write(t.readline())
    out.write(t.readline())
    for line in csv.reader(t,delimiter="\t"):
      gene=line[1].strip()
      pval=genes[gene][1]
      if pval<1:
        o=[str(ranker)]+line[1:]
        o="\t".join(o)+"\n"
        out.write(o)
        ranker+=1

def phevor_in_out(phevor,genes,output):
  out=open(output,'w')
  ranker=0
  with open(phevor) as t:
    out.write(t.readline())
    for line in csv.reader(t,delimiter="\t"):
      gene=line[1].strip()
      pval=genes[gene][1]
      if pval<1:
        o=[str(ranker)]+line[1:]
        o="\t".join(o)+"\n"
        out.write(o)
        ranker+=1

def het_test(info):
  test=False
  if len(info)==1:
    gt=info[0].split(',')[-1]
    if gt=="1":
      test=True
  return test

def main(args):
  genes={}
  exac=tabix.Tabix(args.exac)
  blanked=[1.0,"1.0,1.0","0.0","0.0",""]
  with open(args.pvaast) as t:
    for line in csv.reader(t,delimiter="\t"):
      if line[0]=="RANK":  continue
      if "#" in line[0]:  continue
      line[2]=float(line[2])
      if line[2]==1:
        genes[line[1]]=line[1:]
        continue
      posinfo=line[6:]
      if posinfo==[]:
        modded=[line[1]]+blanked
        genes[line[1]]=modded
        continue
      print line
      print posinfo
      chrom=posinfo[0].split(':')[0]
      print chrom
      gt_freq=genotype_freq(posinfo,exac)
      if args.hets==True and chrom!="chrX":
        het=het_test(posinfo)
        if het==True:
          modded=[line[1]]+blanked
          genes[line[1]]=modded
          continue
      if gt_freq<args.cutoff:
        genes[line[1]]=line[1:]
      else: 
        modded=[line[1]]+blanked
        genes[line[1]]=modded
        continue
  reranked_list=rerank(genes)
  if args.phevor==None and args.grouper==None:
    write_out(args.output,reranked_list)
  elif args.phevor:
    phevor_in_out(args.phevor,genes,args.output)
  elif args.grouper:
    grouper_in_out(args.grouper,genes,args.output)

if __name__=="__main__":
  args=parse_args()
  main(args)
