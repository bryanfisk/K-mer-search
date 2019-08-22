print("importing packages")
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import GridSearchCV
import hashlib
from functions import make_kmer_dict, split_headers, remove_redundant_columns, progress
from scipy import sparse
import numpy as np
import scipy
import time
import os
import sys
import gc

kmer_size = 15
file_dir = "D:/Desktop/metagenomes/reallifesamples"
#file_dir = "D:/Desktop/metagenomes/test_folder"
ref = "D:/Desktop/metagenomes/clostridioides_difficile_genome.fna"

print("loading npz's")
X = sparse.load_npz('X.npz')
y = sparse.load_npz('y.npz')
y = y.toarray()[0]

index_dictionary = {}
with open("index_dictionary.tsv") as file:
    for line in file:
        key, value = list(map(int, line.strip().split()))
        index_dictionary[key] = value

RFC = KNeighborsClassifier(metric = 'l1')
params = {'weights' : ['uniform', 'distance'],'n_neighbors' : [3, 5, 7, 10], 'p' : [2, 4, 6, 8, 10]}
#params = {'weights' : ['uniform']}
GSCV = GridSearchCV(RFC, cv = 5, param_grid = params, n_jobs = 6, verbose = 200)
GSCV.fit(X, y)
print(GSCV.best_params_)
'''
for x in range(10):
    X_train, X_test, y_train, y_test = train_test_split(X, y)
    RFC.fit(X_train, y_train)
    print(RFC.score(X_test, y_test))
'''

#ref_seq = split_headers(ref)
#ref_hash = make_kmer_dict(ref_seq)
#X_val = sparse.lil_matrix((len(os.listdir(file_dir)), 0))
#X_val = sparse.csr_matrix((0, 0))
solution = []
row_count = 0
new_item_count = 1
#columns_in_X_val = sparse.csc_matrix((0, 0))
for file in os.listdir(file_dir):
    pop_keys = []
    input_seq = split_headers(file_dir + '/' + file)
    input_hash = make_kmer_dict(input_seq)
    #max_value = max(index_dictionary.values())
    for key in input_hash.keys():
        if key not in index_dictionary:
            pop_keys.append(key)
            #index_dictionary[key] = max_value + new_item_count
            #new_item_count += 1
    for key in pop_keys:
        input_hash.pop(key)
    #max_key = max(index_dictionary.values())
    row = (0,) * len(input_hash)
    col = tuple(index_dictionary[key] for key in input_hash.keys())
    data = tuple(input_hash.values())
    temp = sparse.csr_matrix((data, (row, col)))
    #if X_val.shape[1] < max_key:
        #X_val.resize(X_val.shape[0], max_key)
        #columns_in_X_val.resize(1, max_key)
    #X_val = sparse.vstack([X_val, temp])
    #columns_in_X_val = columns_in_X_val + temp

#column_list = columns_in_X_val.toarray()[0]


    temp.resize(temp.shape[0], X.shape[1])

    #X_temp = X[:, temp.toarray()[0] > 0]
    #temp = temp[:, temp.toarray()[0] > 0]

    RFC.fit(X, y)
    solution.append(RFC.predict_proba(temp)[0])

with open("output.tsv", "w") as output:
    for sample_index in range(len(os.listdir(file_dir))):
        print("{} : {}".format(os.listdir(file_dir)[sample_index], solution[sample_index]))
        output.write("{}\t{}\n".format(os.listdir(file_dir)[sample_index], solution[sample_index][1]))