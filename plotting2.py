import numpy as np
import sklearn
from sklearn.decomposition import PCA
import random
import matplotlib.pyplot as plt
import pandas as pd
from pandas import DataFrame
from bokeh.io import output_notebook, show
from bokeh.plotting import figure
from bokeh.models import TapTool, CustomJS, ColumnDataSource
from bokeh.models import HoverTool
from bokeh.io import output_file
from bokeh import colors


f = open('3popreal.eigenvec')
q=open('superlabels.txt')

superlabels={}
for line in q:
    line=line.strip().split(" ")
    superlabels[line[0]] = [ line[1],line[2]]


   
data=[]
labels=[]
samplpelabels=[]
for line in f:
    line=line.strip().split(' ')
    for j in range(2,len(line)):
        line[j]=float(line[j])
    data.append(line[2:])
    labels.append(line[0])
    samplpelabels.append(line[1])

mysuperlabels=[ " " for x in labels]
mysupersuperlabels=[ " " for x in labels]
for r in range(0,len(labels)):
    mysuperlabels[r]=superlabels[labels[r]][1]
    mysupersuperlabels[r]=superlabels[labels[r]][0]


output_file('plot.html')
source = ColumnDataSource(
        data=dict(
            x=[x[0] for x in data],
            y=[x[1] for x in data],
            label=labels,
            ID=samplpelabels,
            superlabel=mysuperlabels,
            supersuperlabel=mysupersuperlabels,
            colour=[ colors.RGB( 255/(len(set(mysupersuperlabels)))*(list(set(mysupersuperlabels)).index(mysupersuperlabels[x])),255/(len(set(mysuperlabels)))*(list(set(mysuperlabels)).index(mysuperlabels[x])) ,255/(len(set(labels)))*(list(set(labels)).index(labels[x])) ) for x in range(0,len(labels))]
        )
    )


hover = HoverTool(
        tooltips=[
            ("ID", "@ID"),
            ("population", "@label"),
            ("Area", "@superlabel"),
            ("Continent", "@supersuperlabel")
        ]
    )

p = figure(plot_width=1300, plot_height=700, tools=[hover], title="Pvalue of sudies through time")

p.circle('x', 'y', source=source,fill_color="colour")

show(p)
