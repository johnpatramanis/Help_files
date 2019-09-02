#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 13:52:01 2018

@author: aggelos
"""

#general
import math
import numpy as np
from numpy import genfromtxt
import timeit
import matplotlib.pyplot as plt
import re
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

#prediction
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report
from sklearn.metrics import accuracy_score
from sklearn.metrics import make_scorer
from sklearn.metrics import roc_curve


#feature selection
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import SelectFromModel
from sklearn.feature_selection import mutual_info_classif

#how to crash a server 101
from sklearn.externals import joblib
from sklearn.externals.joblib import Parallel, delayed
import multiprocessing


start_time = timeit.default_timer()

nested_scores = {}

clf_names = ['neutrality', 'selection']

true_pos = []
false_pos = []

accuracies = np.zeros((61,3)).astype(float)
print(accuracies.shape)

dtst = 0
def get_the_rep(y_true, y_pred):
    global dtst
    nested_scores.update({'params'+str(dtst) : classification_report(y_pred,y_true,target_names = clf_names)})
    
    fpr, tpr, _ = roc_curve(y_true, y_pred, pos_label=0)
    true_pos.append(tpr)
    false_pos.append(fpr)         
    
    
    dtst = dtst + 1
    return accuracy_score(y_true, y_pred)


#def get_the_rep(y_true, y_pred):
#    global dtst
#    best_params.update({'params'+str(dtst) : classification_report(y_pred,y_true,target_names = clf_names)})
#    dtst = dtst + 1
#    return accuracy_score(y_true, y_pred)



#------------------------- Data Manipulation -------------------------
#seed = np.array([1])
#seed.rand(1)
#print (seed)
path = 'summ_stats/'

chosen_sets = [7, 8, 9, 10, 11, 12, 19, 20, 21, 22, 23, 24, 30, 31,
               32, 33, 35, 36, 43, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 60]

neut_test = {}
slct_test = {}

neut_data = {}
slct_data = {}


num_sets  = 61 #number of datasets + 1

for i in range (1,num_sets):
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


#------------------ Support Vector Machines ---------------
print("Support Vector Machines")
pipeline_svm = make_pipeline(StandardScaler(), SVC())
cst = [1, 2, 5, 7, 10] #starting from hard-margin and loosening to see performance
ker = ['poly']
deg = [1] #just for the polynomial kernel


hyperparameters = { 'svc__C' : cst,
                   'svc__kernel': ker,
                   'svc__degree': deg,
                   }

inner_cv = KFold(n_splits=10, shuffle=True)
outer_cv = KFold(n_splits=5, shuffle=True)

## Pass the gridSearch estimator to cross_val_score
clf = GridSearchCV(pipeline_svm, param_grid=hyperparameters, cv=inner_cv, n_jobs=-1)



#for each demographic pair we perform nested cross validation
for feats in range (40, 41):
    
    set_scores_svm = []
 #   global true_pos, false_pos, nested_scores

    nested_scores = {}
    true_pos = []
    false_pos = []
    
    print (feats)
    for i in range(1, num_sets):
        print(i)
        data = np.concatenate((neut_data['neut'+str(i)],slct_data['sel'+str(i)]), axis = 0)
        selected = SelectKBest(score_func=mutual_info_classif, k=feats).fit(data, labels[:,0])
        current = selected.transform(data)
        
        kappa = cross_val_score(clf, X=current, y=labels[:,0], cv=outer_cv, scoring = make_scorer(get_the_rep)).mean()
        set_scores_svm.append(kappa)    
    
    with open ("split_res/acc/tot_acc_split" + str(feats) + ".txt", "w") as rf:
        for j in range (0, 60):
           rf.write(str(set_scores_svm[j]) + "\n")
    
    np.savetxt("split_res/tpr/tpr_split" + str(feats) + ".txt", np.array(true_pos[:]),  fmt = '%1.8f')
    np.savetxt("split_res/fpr/fpr_split" + str(feats) + ".txt",  np.array(false_pos[:]),  fmt = '%1.8f' )
    #np.savetxt("split_res/nested/rep_split" + str(feats) + ".txt", np.array(true_pos[:]),  fmt = '%1.8f')
    	#rp.write(str(nested_scores))
