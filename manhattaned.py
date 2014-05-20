import matplotlib.pyplot as plt
import numpy as np
import csv
import random

chroms=list(np.arange(1,23))+['X','Y']
chroms=[str(i) for i in chroms]

def rotate(nist):
	return nist[1:]+nist[:1]

scores={}
with open("WPW_Fam1_pvaast_indel_3.phevor.txt") as t:
	for line in csv.reader(t,delimiter='\t'):
		if "#" in line[0]:
			continue
		else:
			scores[line[1]]=float(line[2])

genes={}
last=[0]
with open("manhattanSorted.bed") as q:
	for line in csv.reader(q,delimiter='\t'):
		glist=line[3].split(',')
		for i in glist:
			if line[0] not in genes.keys():
				genes[line[0]]={}
			genes[line[0]][i]={}
			if i in scores.keys():
				genes[line[0]][i]['score']=scores[i]
			else:
				genes[line[0]][i]['score']=0
			coord=random.randint(int(line[1]),int(line[2]))
			genes[line[0]][i]['coord']=coord


for i in range(len(chroms)):
	alcor=[]
	for j in genes[chroms[i]].keys():
		alcor.append(genes[chroms[i]][j]['coord'])
	if len(last) > 1:
		last.append(max(alcor)+last[i])
	else:
		last.append(max(alcor))

corr=[]
scor=[]
t=0
sig=['MYH6','LAMB2','FUZ']
sigd={}
for k in chroms:
	for j in genes[k].keys():
		if j in sig:
			sigd[j]=[genes[k][j]['coord']+last[t], genes[k][j]['score']]
		corr.append(genes[k][j]['coord']+last[t])
		scor.append(genes[k][j]['score'])
	t+=1

print sigd
centering=[(last[i]+last[i+1])/2 for i in range(len(last)) if i < len(last)-1]
colors=['b','g','r','y']

plt.figure(figsize=(16,16))
ax=plt.subplot()
ax.set_ylim([-0.05,6])
ax.set_xlim([min(corr),max(corr)])
ax.set_xticks(centering)
ax.set_xticklabels(chroms)
plt.xlabel("Chromosome",fontsize=16)
plt.ylabel("Phevor Score",fontsize=16)
plt.title("Candidate Genes Re-Ranked by Phevor",fontsize=24)
t=0
for k in chroms:
	corr=[]
	scor=[]
	for j in genes[k].keys():
		corr.append(genes[k][j]['coord']+last[t])
		scor.append(genes[k][j]['score'])
	plt.scatter(corr,scor,c=colors[0],marker='.',s=500,edgecolors='none')
	colors=rotate(colors)
	t+=1
plt.text(sigd['MYH6'][0]+20000000,sigd['MYH6'][1],"MYH6",fontsize=16)
plt.show()
