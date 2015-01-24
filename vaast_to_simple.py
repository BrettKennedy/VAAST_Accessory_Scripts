#!/usr/bin/env python

import csv, sys
import operator
csv.field_size_limit(sys.maxsize)
genes={}
lod=0
pval=0
gene=0
score=0
ci=0

out=open("trying.simple",'w')
with open("FamAF_Fam3_pvaast_snv_2.vaast") as t:
    for line in csv.reader(t,delimiter="\t"):
        if ">" in line[0]:
            new_gene=line[1]
        if "TU" in line[0] or "TR" in line[0]:
            tsco=float(line[1].split('(')[0].strip())
            if tsco>0:
                score=tsco
        if "LOD_SCORE" in line[0]:
            lod=line[0].split(':')[1].split(',')[0]
        if "genome_permutation_p" in line[0]:
            pval=float(line[0].split(':')[1])
        if "genome_permutation_0.95_ci" in line[0]:
            ci=line[0].split(':')[1]
        if "RANK" in line[0]:
            simple_info=[new_gene,pval,ci,score,lod,"blank"]
            if new_gene not in genes:
                genes[new_gene]=simple_info
            elif pval<genes[new_gene][1]:
                genes[new_gene]=simple_info

list_to_sort=[]
for g in genes:
    list_to_sort.append(genes[g])
    
    
lis=sorted(list_to_sort, key=operator.itemgetter(1))


out.write("RANK\tGene\tp-value\tp-value-ci\tScore\tLOD\tVariants\n")
rank=1
for g in lis:
    g=[rank]+g
    info=[str(i) for i in g]
    info="\t".join(info)+"\n"
    out.write(info)
    rank+=1
out.close()
