#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import csv, argparse, random

__author__="Brett J. Kennedy and Gordon Lemmon"
__date__="January 24 2015"

plt.rcParams['figure.autolayout'] = True

def get_gene_list(path):
	with open(path) as fh:
		genes=[]
		for line in fh:
			genes.append( line.strip() )
	return genes

from os.path import join

gene_type_dir='../gene_categories'
genes=dict(
		cilium=get_gene_list(join(gene_type_dir, 'Cilium_Genes.txt')),
		chrX=get_gene_list(join(gene_type_dir,'chrX_Genes.txt')),
		chromatin=get_gene_list(join(gene_type_dir,'chromatin_genes.txt'))
)

def add_gene_categories(info):
	categories = []
	targets = info['targets']
	for category in ['cilium','chrX','chromatin']:
		for target in targets:
			if target in genes[category]: categories.append(category)
	info['categories']=categories	
	return info

def get_colors_from_category(category):
	colors=[]
	if category == 'cilium': colors.append('LightSkyBlue')
	if category == 'chrX': colors.append('LightSalmon')
	if category == 'chromatin': colors.append('LightGreen')
	return colors

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

def parse_scores(path,phevor):
	"""Parses either .simple or .phevor output and returns dictonary\
	of sores by gene."""
	scores={}
	highs=0
	high_gene=None

	start=False
	for line in open(path):
		line=line.strip()
		if len(line)==0: continue
		if not start:
			if line.startswith('# RANK') or line.startswith('RANK'):
				start=True
			continue
		if line.startswith('#'): continue
		line=line.split()
		gene=line[1]
		HPO=0.0 # for pvaast
		score=float(line[2])
		if score < 0: score = 0
		if score > 0:
			if score > highs:
				highs=score
				high_gene=gene
			if phevor:
				terms=line[5].split('|')
				try: HPO=float(terms[1]) # should crash on parse error
				except: pass
			else: #if not path.endswith('.renorm'):
				score=-np.log10(score)
		scores[gene]=(score,HPO) 
	return scores,highs,high_gene

def populate_scores_coords(gff3,scores):
    """Parses the gff3 and creates a dictionary of genomic\
     coordinates 'coord' and 'scores'"""
    genes={}
    with open(gff3) as q:
        for line in csv.reader(q,delimiter="\t"):
            chrom=line[0]
            if "#" in chrom:  continue
            if line[2]!="gene": continue
            info=parse_info(line[8])
            name=info['Name']
            start=int(line[3])
            stop=int(line[4])
            if chrom not in genes.keys():
                genes[chrom]={}
            genes[chrom][name]={}
            if name in scores.keys():
                genes[chrom][name]['score']=scores[name][0]
                genes[chrom][name]['color']=scores[name][1]
            else:
                genes[chrom][name]['score']=0.0
                genes[chrom][name]['color']=0.0
            coord=random.randint(start,stop)
            genes[chrom][name]['coord']=coord
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

def add_colors(ax, categories):
	colors=[]
	for category in categories:
		colors.extend(get_colors_from_category(category))
	if not len(colors): return
	XMIN,XMAX = ax.get_xlim()
	YMAX = ax.get_ylim()[1]
	xsize = (XMAX - XMIN)/len(colors)
	for i,color in enumerate(colors):
		xmin = XMIN + i*xsize
		xmax = xmin + xsize
		ax.bar(xmin,YMAX,xsize,color=color,zorder=1,linewidth=0)

def panel_key_fxn(info):
	types=('cilium','chrX','chromatin')
	types_dict=dict(zip(types,[1,2,3]))
	categories = info['categories']
	if len(categories) == 0:
		return 0
	if len(categories) == 3:
		return 6

	return np.mean([types_dict[category] for category in categories]) 

def draw_manhattans(
		inputs,
		output,
		gff3,
		width=15,
		height=10,
		point_size=300,
		columns=1,
		rows=1,
		phevor=True,
):
	chroms=list(np.arange(1,23))+['X','Y']
	chroms=["chr"+str(i) for i in chroms]

	plt.figure(figsize=(width,height))

	i = 0
	highest=0

	inputs = [add_gene_categories(info) for info in inputs]

	inputs.sort(key=lambda x: int(x['name']))
	#inputs.sort(key=panel_key_fxn)
	for info in inputs:
		scores,highS,high_gene=parse_scores(info['path'],phevor)
		if highS > highest:
			highest=highS
		info['highS']=highS
		info['high_gene']=high_gene
		info['genome']=populate_scores_coords(gff3,scores)

	for info in inputs:	
		i += 1
		last=last_list(info['genome'],chroms)
		ax=plt.subplot(rows, columns, i)
		ax.set_ylim(bottom=0, top=highest*1.1)
		ax.set_xlim(0,max(last) )
		plt.setp(ax.get_xticklabels(), visible=False)
		add_colors(ax, info['categories'])
		draw_manhattan(
				info['highS'],
				info['high_gene'],
				info['genome'],
				info['targets'],
				last,
				None,
				None,
				' '.join([info['name']]+list(set(info['diseases']))),
				width,
				height,
				point_size,
				columns,
				rows,
				chroms,
		)
		#ax.text(500,y_max*1.1*0.80,info['name'],fontsize=24,color='black',weight='bold')
		ax.axes.get_xaxis().set_visible(False)
	plt.savefig(output,dpi=300)
	plt.clf()

def draw_manhattan(
		highS,
		high_gene,
		genome,
		targets,
		last,
		xlab="Position",
		ylab="Score",
		title="pVaast or Phevor output",
		width=15,
		height=5,
		point_size=300,
		columns=1,
		rows=1,
		chroms=[],
		plot_all_above_zero = True
):

	corr=[]
	scor=[]

	centering=[(last[i]+last[i+1])/2 for i in range(len(last)) if i < len(last)-1]

	#colors=['k'] # black only 
	colors=['b','g','r','y'] 

	if xlab: plt.xlabel(xlab,fontsize=12)
	if ylab: plt.ylabel(ylab,fontsize=16)
	if title:	plt.title(title,fontsize=24)
	t=0
	for c in chroms:
		zero_coor=[]
		zero_scor=[]
		zero_colors=[]
		corr=[]
		scor=[]
		textX=[]
		target_text=None
		
		chrom = genome[c]
		for gene in chrom:
			if chrom[gene]['score'] > 0:
				corr.append(chrom[gene]['coord']+last[t])
				scor.append(chrom[gene]['score'])
				if gene in targets:
					target_text=(gene,chrom[gene]['coord']+last[t],chrom[gene]['score'])
				else:
					textX.append((gene,chrom[gene]['coord']+last[t],chrom[gene]['score']))
			else:
				if gene in targets:
					target_text=(gene,chrom[gene]['coord']+last[t],chrom[gene]['score'])
				zero_coor.append(chrom[gene]['coord']+last[t])
				zero_scor.append(chrom[gene]['score']) # shouldn't this just be zero?
		plt.scatter(zero_coor,zero_scor,marker='.',s=300,color=colors[0],edgecolors='none',zorder=2)
		plt.scatter(corr,scor,marker='.',s=point_size,color=colors[0],edgecolors='black',zorder=3)

		for text in textX:
			gene,x,y = text
			if gene == high_gene:
				plt.text(x,y,gene,fontsize=13,color=colors[0],horizontalalignment='center')
			if plot_all_above_zero:
				plt.text(x,y,gene,fontsize=13,color=colors[0],horizontalalignment='center')
		if target_text:
			plt.text(target_text[1],target_text[2],target_text[0],color='k',horizontalalignment='center',weight='bold')

		colors=rotate(colors)
		t+=1
