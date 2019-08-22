print("importing packages")
import argparse
from functions import dir_to_matrix
from scipy import sparse
import csv
from sklearn.utils import resample
import os
from sklearn.neighbors import KNeighborsClassifier

file_dir_train = "D:/Desktop/metagenomes/real_life_train/train"
file_dir_test = "D:/Desktop/metagenomes/real_life_train/test"
reference = "D:/Desktop/metagenomes/clostridioides_difficile_genome.fna"

X_train, y_train, index_dictionary_train = dir_to_matrix(file_dir_train, reference)
X_test, y_test, index_dictionary_test = dir_to_matrix(file_dir_test, reference)
sample_names_train = os.listdir(file_dir_train)
sample_names_test = os.listdir(file_dir_test)

X_train, y_train, sample_names_train = resample(X_train, y_train, sample_names_train, n_samples = 100)
X_test, y_test, sample_names_test = resample(X_test, y_test, sample_names_test, n_samples = 100)

print(X_train.shape, y_train.shape, sample_names_train.shape)
print(X_test.shape, y_test.shape, sample_names_test.shape)

KNN = KNeighborsClassifier()
KNN.fit(X_train, y_train)
print(KNN.score(X_test, y_test))