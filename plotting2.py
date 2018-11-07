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

data=[]
familylabels=[]
idlabels=[]
for line in f:
    line=line.strip().split()
    familylabels.append(line[0])
    idlabels.append(line[1])
    datalabels.append([line[2],line[3]])



output_file('plot.html')
source = ColumnDataSource(
        data=dict(
            x=[x[0] for x in data],
            y=[x[1] for x in data],
            label=idlabels,
            family=familylabels,
            colour=['blue'*len(x)]
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
