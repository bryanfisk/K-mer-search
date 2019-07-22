print("importing packages")
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
import hashlib
from functions import make_kmer_dict, split_headers, remove_redundant_columns, progress
from scipy import sparse
import numpy as np
import scipy
import time
import os

kmer_size = 15
#file_dir = "D:/Desktop/metagenomes/reallifesamples"
file_dir = "D:/Desktop/metagenomes/test_folder"
ref = "D:/Desktop/metagenomes/clostridioides_difficile_genome.fna"

print("loading npz's")
X = sparse.load_npz('X_train.npz')
y = sparse.load_npz('y_train.npz')
y = y.toarray()[0]

index_dictionary = {}
with open("index_dictionary.tsv") as file:
    for line in file:
        key, value = list(map(int, line.strip().split()))
        index_dictionary[key] = value

ref_seq = split_headers(ref)
ref_hash = make_kmer_dict(ref_seq)
#X_val = sparse.lil_matrix((len(os.listdir(file_dir)), 0))
X_val = sparse.csr_matrix((0, 0))
row_count = 0
new_item_count = 1
columns_in_X_val = sparse.csc_matrix((0, 0))
for file in os.listdir(file_dir):
    input_seq = split_headers(file_dir + '/' + file)
    input_hash = make_kmer_dict(input_seq)
    max_value = max(index_dictionary.values())
    input_hash_length = len(input_hash)
    for key in input_hash.keys():
        if key not in index_dictionary:
            index_dictionary[key] = max_value + new_item_count
            new_item_count += 1
    max_key = max(index_dictionary.values())
    row = [0 for x in range(len(input_hash))]
    col = [index_dictionary[key] - 1 for key in input_hash.keys()]
    data = np.fromiter(input_hash.values(), dtype = int) # running out of memory, wat do?
    temp = sparse.csr_matrix((data, (row, col)))
    if X_val.shape[1] < max_key:
        X_val.resize(X_val.shape[0], max_key)
        columns_in_X_val.resize(1, max_key)
    X_val = sparse.vstack([X_val, temp])
    columns_in_X_val = columns_in_X_val + temp

column_list = columns_in_X_val.toarray()[0]

'''
column_list = []
for column in range(columns_in_X_val.shape[1]):
    if columns_in_X_val[:, column] > 0:
        column_list.append(column)
'''

if X.shape[1] < X_val.shape[1]:
    X.resize(X.shape[0], X_val.shape[1])
elif X.shape[1] > X_val.shape[1]:
    X_val.resize(X_val.shape[0], X.shape[1])

#keepset = remove_redundant_columns(X_val)
X = X[:, column_list > 0]
X_val = X_val[:, column_list > 0]

X_train, X_test, y_train, y_test = train_test_split(X, y)

RFC = RandomForestClassifier(n_estimators = 100)
RFC.fit(X_train, y_train)
solution = RFC.predict_proba(X_val)

with open("output.csv", "w") as output:
    for sample_index in range(len(os.listdir(file_dir))):
        print("{} : {}".format(os.listdir(file_dir)[sample_index], solution[sample_index]))
        output.write("{},{}\n".format(os.listdir(file_dir)[sample_index], solution[sample_index][1]))






'''
#split dictionary into list of sub-dictionaries
def split_dict(dict, num_splits):
dict_list = []
for split in range(num_splits):
    split_size = math.ceil(len(dict)/num_splits)
    temp_dict = dict(list(dict.items())[split * split_size : split * split_size + split_size])
    dict_list.append(temp_dict)
while {} in dict_list:
    dict_list.remove({})
return dict_list
'''