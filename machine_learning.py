print("importing packages")
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import MultinomialNB
import hashlib
from scipy import sparse
import numpy as np
import time
import sys
import os
import re

kmer_size = 15
real_samples = "D:/Desktop/metagenomes/reallifesamples"
starttime = time.time()
reference = "D:/Desktop/metagenomes/clostridioides_difficile_genome.fna"

print("loading npz's")
X = sparse.load_npz('X_train.npz')
y = sparse.load_npz('y_train.npz')
y = y.toarray()[0]

X_train, X_test, y_train, y_test = train_test_split(X, y)
MNB = MultinomialNB()
MNB.fit(X_train, y_train)
print(MNB.score(X_test, y_test))