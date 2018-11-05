import numpy as np
import matplotlib.pyplot as plot
import colorsys
import getopt
import sys, pdb
import re
import argparse

parser = argparse.ArgumentParser()  
parser.add_argument('--file',nargs=1,type=str)
parser.add_argument('--labels',nargs=1,type=str)

args = parser.parse_args()





def plot_admixture(admixture, population_indices, population_labels, title):

    N,K = admixture.shape
    colors = [colorsys.hsv_to_rgb(h,0.9,0.7) for h in np.linspace(0,1,K+1)[:-1]]
    text_color = 'k'
    bg_color = 'w'
    fontsize = 12

    figure = plot.figure(figsize=(5,3))

    xmin = 0.13
    ymin = 0.2
    height = 0.6
    width = 0.74
    indiv_width = width/N
    subplot = figure.add_axes([xmin,ymin,width,height])
    [spine.set_linewidth(0.001) for spine in subplot.spines.values()]

    for k in range(K):
        if k:
            bottoms = admixture[:,:k].sum(1)
        else:
            bottoms = np.zeros((N,),dtype=float)

        lefts = np.arange(N)*indiv_width
        subplot.bar(lefts, admixture[:,k], width=indiv_width, bottom=bottoms, facecolor=colors[k], edgecolor=colors[k], linewidth=0.4)

        subplot.axis([0, N*indiv_width, 0, 1])
        subplot.tick_params(axis='both', top=False, right=False, left=False, bottom=False)
        xtick_labels = tuple(map(str,['']*N))
        subplot.set_xticklabels(xtick_labels)
        ytick_labels = tuple(map(str,['']*K))
        subplot.set_yticklabels(ytick_labels)

    position = subplot.get_position()
    title_position = (0.5, 0.9)
    figure.text(title_position[0], title_position[1], title, fontsize=fontsize, \
        color='k', horizontalalignment='center', verticalalignment='center')

    for p,popname in enumerate(population_labels):
        indices = np.where(population_indices==p)[0]
        if indices.size>0:
            vline_pos = (indices.max()+1)*indiv_width 
            subplot.axvline(vline_pos, linestyle='-', linewidth=0.2, c='#888888')
            label_position = (xmin+(2*indices.min()+indices.size)*0.5*indiv_width, ymin-0.01)
            figure.text(label_position[0], label_position[1], popname, fontsize=6, color='k', \
                horizontalalignment='right', verticalalignment='top', rotation=70)

    return figure

	
file=open(str(args.file[0]))
admixturedata=[]
for l in file:
	datapoint=l.strip().split(" ")
	admixturedata.append(datapoint)
admixturedata=np.asarray(admixturedata)

labels=open(str(args.labels[0]))
mylabels=[]
for line in labels:
	labelpoint=line.strip().split(" ")
	mylabels.append(labelpoint[0])

popindices=[]
counter=0
for i in range(1,len(mylabels)):
	counter+=1
	if mylabels[i]!=mylabels[i-1]:
		popindices.append(counter)

plot_admixture(admixturedata,popindices,mylabels,"ADMIXTURE ANALYSIS")
plot.show()