#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import csv, argparse, random

def rotate(nist):
    """Spins a list"""
    return nist[1:]+nist[:1]

def parse_info(info):
    """creates a dictionary of the info field in a GFF3"""
    infoD={}
    info=info.split(';')
    for e in info:
            e=e.split('=')
            try:    infoD[e[0]]=e[1]
            except: continue
    return infoD

def parse_scores(input,phevor2):
    """Parses either .simple or .phevor output and returns dictonary\
    of sores by gene."""
    scores={}
    highs=0
    if phevor2==False:
        with open(input) as t:
            for line in csv.reader(t,delimiter="\t"):
                if line[0]=="RANK": continue
                sco=-np.log10(float(line[2]))
                scores[line[1].strip()]=sco
                if sco>highs:
                    highs=sco
    
    if phevor2==True:
        with open(input) as t:
            for line in csv.reader(t,delimiter="\t"):
                if "#" in line[0]:  continue
                sco=float(line[2])
                scores[line[1].strip()]=sco
                if sco>highs:
                    highs=sco
    return scores,highs

def populate_scores_coords(gff3,scores):
    """Parses the gff3 and creates a dictionary of genomic\
     coordinates 'coord' and 'scores'"""
    genes={}
    with open(gff3) as q:
        for line in csv.reader(q,delimiter="\t"):
            if "#" in line[0]:  continue
            if line[2]!="gene": continue
            info=parse_info(line[8])
            name=info['Name']
            chrom=line[0]
            start=int(line[3])
            stop=int(line[4])
            if line[0] not in genes.keys():
                genes[line[0]]={}
            genes[line[0]][name]={}
            if name in scores.keys():
                genes[line[0]][name]['score']=scores[name]
            else:
                genes[line[0]][name]['score']=0
            coord=random.randint(start,stop)
            genes[line[0]][name]['coord']=coord
    return genes

def last_list(genes,chroms):
    last=[0]
    for i in range(len(chroms)):
    	alcor=[]
    	for j in genes[chroms[i]].keys():
    		alcor.append(genes[chroms[i]][j]['coord'])
    	if len(last) > 1:
    		last.append(max(alcor)+last[i])
    	else:
    		last.append(max(alcor))
    return last

def parse_args():
    parser=argparse.ArgumentParser(description="Creates a manhattan plot from VAAST or Phevor output")
    parser.add_argument("input", help="simple vaast or pVAAST output")
    parser.add_argument("output", help="save file to this handle",type=str)
    parser.add_argument("gff3", help="bed file of genomic coordinates of genes")
    parser.add_argument("--genes",help="comma seperated genes of interest to label in plot",default=None)
    parser.add_argument("--phevor2",help="use this arguemnt if the input is phevor2 rather than VAAST/pVAAST",action='store_true',default=False)
    return parser.parse_args()

# genes is a gene name or list of names
def plot_stuff(input,output,gff3,sig,phevor2):
    if not isinstance(sig, list):
      sig = [sig]
    chroms=list(np.arange(1,23))+['X','Y']
    chroms=["chr"+str(i) for i in chroms]
    
    scores,highS=parse_scores(input,phevor2)
    genes=populate_scores_coords(gff3,scores)
    last=last_list(genes,chroms)
    
    corr=[]
    scor=[]
    t=0
    sigd={}
    for k in chroms:
    	for j in genes[k].keys():
    		if j in sig:
    			sigd[j]=[genes[k][j]['coord']+last[t], genes[k][j]['score']]
    		corr.append(genes[k][j]['coord']+last[t])
    		scor.append(genes[k][j]['score'])
    	t+=1

    #print sigd
    centering=[(last[i]+last[i+1])/2 for i in range(len(last)) if i < len(last)-1]
    colors=['b','g','r','y']

    plt.figure(figsize=(15,10))
    ax=plt.subplot()
    ax.set_ylim([-0.05,highS+1])
    ax.set_xlim([min(corr),max(corr)])
    ax.set_xticks(centering)
    ax.set_xticklabels(chroms,rotation=45)
    plt.xlabel("Chromosome",fontsize=16)
    if phevor2==True:
        ylab="Phevor Score"
        title="Genes Re-Ranked by Phevor"
    else:
        ylab="$-log_{10}$ pVAAST p-value"
        title="Genes Scored by pVAAST"
    plt.ylabel(ylab,fontsize=16)
    plt.title(title,fontsize=24)
    t=0
    for k in chroms:
    	corr=[]
    	scor=[]
    	for j in genes[k].keys():
    		corr.append(genes[k][j]['coord']+last[t])
    		scor.append(genes[k][j]['score'])
    	plt.scatter(corr,scor,c=colors[0],marker='.',s=300,edgecolors='none')
    	colors=rotate(colors)
    	t+=1
    if args.genes!=None:
        for g in sig:
            plt.text(sigd[g][0]+20000000,sigd[g][1],g,fontsize=16)
    plt.savefig(args.output,dpi=300)
    plt.show()

if __name__=="__main__":
    args = parse_args()
    genes=[]
    if args.genes!=None:
        genes=args.genes.split(',')
    plot_stuff(args.input,args.output,args.gff3,genes,args.phevor2)
    
