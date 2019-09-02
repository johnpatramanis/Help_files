# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 14:58:03 2018

@author: aaggelos
"""

import math
import numpy as np
from numpy import genfromtxt
from numpy import random
import timeit
import matplotlib.pyplot as plt
#import pandas as pd #not sure if needed
#model selection
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
#classifiers
from sklearn.naive_bayes import MultinomialNB
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
#cross validation
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import KFold
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_score
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.externals import joblib

from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score

start_time = timeit.default_timer()




#------------------------- Data Manipulation -------------------------
#seed = np.array([1])
#seed.rand(1)
#print (seed)
path = 'summ_stats/'
neut_test = {}
slct_test = {}

neut_data = {}
slct_data = {}

labels = np.zeros((1000,1)).astype(int)*(-1)
for i in range (1,61):
    if (i < 10):
        neutpath = path + 'neutral0' + str(i) + '.txt'
        selpath = path + 'sel0' + str(i) + '.txt'
    else:
        neutpath = path + 'neutral' + str(i) + '.txt'
        selpath = path + 'sel' + str(i) + '.txt'
    neut = genfromtxt(neutpath, delimiter='\t', dtype='float')
    sel = genfromtxt(selpath, delimiter = '\t', dtype='float')
    #X_train, X_test, Y_train, Y_test = train_test_split(neut_data, labels, test_size=0.2,random_state=123,stratify=labels)
    neut = np.delete(neut, (0), axis=0)
    sel = np.delete(sel, (0), axis=0)
    neut_data.update({'neut'+str(i) : neut})
    slct_data.update({'sel'+str(i) : sel})
    
print(neutpath)
print(selpath)

print(neut_data['neut1'])

smp = 1000 #number of samples for each class
neutral_labels = np.ones((1,smp)).astype(int)
print(neutral_labels.shape)
selection_labels = np.zeros((1,smp)).astype(int)*(-1)
print(selection_labels.shape)
labels = np.concatenate((neutral_labels, selection_labels), axis = 1).T
c,r = labels.shape
labels.reshape(c,)
print(labels.shape)

#-------------------------- Random Forests --------------------------
print("Random Forest")


pipeline_rf = make_pipeline(StandardScaler(), RandomForestClassifier( n_jobs=-1))

crit = ['gini', 'entropy']
numt = [100, 150] #number of trees
feats = ['log2','sqrt'] #max_features considered during each split
depth = [ 20, 30, None] #maximum depth each tree is allowed to reach
# 
hyperparameters = {#'randomforestclassifier__criterion' : crit,
                   'randomforestclassifier__n_estimators' : numt,
                   'randomforestclassifier__max_features' : feats,
                   'randomforestclassifier__max_depth': depth}


nested_score = np.zeros([60, 1]).astype(int)
print(nested_score.shape)

inner_cv = KFold(n_splits=10, shuffle=True)
outer_cv = KFold(n_splits=5, shuffle=True)

## Pass the gridSearch estimator to cross_val_score
clf = GridSearchCV(pipeline_rf, param_grid=hyperparameters, cv=inner_cv)

set_scores = []

#for each demographic pair we perform nested cross validation

for i in range(1, 61):
    data = np.concatenate((neut_data['neut'+str(i)],slct_data['sel'+str(i)]), axis = 0)
    kappa = cross_val_score(clf, X=data, y=labels[:,0], cv=outer_cv, n_jobs=-1).mean()
    nested_score[i-1,0] = kappa
    set_scores.append(kappa)
    print("dataset" + str(i))
with open ("forest_res.txt",'w') as rf:
    for j in range (0, 60):
        rf.write(str(set_scores[j]))
    
print(kappa)

print(set_scores)
