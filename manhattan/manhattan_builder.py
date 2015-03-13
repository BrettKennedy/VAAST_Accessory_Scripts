#!/usr/bin/env python

from glob import glob
from os.path import join, basename
from manhattaned_grid import draw_manhattans

out='results/manhattan/png/brett_results'
gff3='../../data/ref_GRCh37.p13_14-09-20.genes_only.gff3'

def parse_targets(target_file='inputs/Phevor.info.noheader'):
	targets = dict()
	diseases = dict()
	for line in open( target_file ):
		ped,disease,gene = line.strip().split()
		if ped in targets:
			targets[ped].append(gene)
			diseases[ped].append(disease)
		else:
			targets[ped]=[gene]
			diseases[ped]=[disease]
	return targets,diseases

targets,diseases=parse_targets()

def process_glob(glob_me, png, phevor):
		peds=[]
		for path in glob(glob_me):
			base = basename(path)
			name = base.split('.')[0]
			peds.append({'name':name,'path':path,'targets':targets[name],'diseases':diseases[name] })
		#if input_type == 'pvaast':
		#	phevor=False
		draw_manhattans(peds,png,gff3,width=15,height=20,point_size=300,columns=6,rows=10, phevor=phevor)

def main():
	count=0
	phevor_dir = join('results','from_brett')
	png = join(out,'grouper.png')
	process_glob(join(phevor_dir,'grouper','*.grouper'), png, True)
	png = join(out,'grouper.exac.png')
	process_glob(join(phevor_dir,'grouper_filtered','*.grouper.exac'), png, True)
	png = join(out,'phevor_renorm.png')
	process_glob(join(phevor_dir,'phevor_norm','*.phevor2.renorm'), png, True)

if __name__ == '__main__':
	main()
