import numpy as np
import argparse
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


parser = argparse.ArgumentParser()  
parser.add_argument('--eig',nargs='+',type=str)
parser.add_argument('--out',nargs='+',type=str)




args = parser.parse_args()

f = open(args.eig[0])




data=[]
familylabels=[]
idlabels=[]
for line in f:
    line=line.strip().split()
    familylabels.append(line[0])
    idlabels.append(line[1])
    data.append([line[2],line[3]])

uniqlabels=list(set(familylabels))
colordic={}
for k in uniqlabels:
    colordic[k]=colors.RGB(random.randint(1,255),random.randint(1,255),random.randint(1,255))
print(colordic) 
colorz=[colordic[x] for x in familylabels]  
    
    


output_file('{}.html'.format(args.out[0]))
source = ColumnDataSource(
        data=dict(
            x=[float(x[0]) for x in data],
            y=[float(x[1]) for x in data],
            label=idlabels,
            family=familylabels,
            colour=colorz
        )
    )


hover = HoverTool(
        tooltips=[
            ("ID", "@label"),
            ("population", "@family")
        ]
    )

p = figure(plot_width=1300, plot_height=700, tools=[hover], title="Population PCA")

p.circle('x', 'y', source=source,fill_color="colour")

show(p)
