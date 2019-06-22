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
X = sparse.load_npz('extracted_features.npz')
y = np.load('y_values.npz')
y = y['arr_0']

X_train, X_test, y_train, y_test = train_test_split(X, y)
MNB = MultinomialNB()
MNB.fit(X_train, y_train)
print(MNB.score(X_test, y_test))

#progress bar
def progress(count, total, status = ''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))
    percents = round( 100.0 * count / float(total), 1)
    bar = '■' * filled_len + '-' * (bar_len - filled_len)
    sys.stdout.write('|{:}| {:}{:}  {:}\r'.format(bar, percents, '%', status))
    sys.stdout.flush()

#split headers from sequences from fastq file
def split_headers(file):
    pattern = re.compile('.*\.fastq$')
    headers = []
    sequences = []
    with open(file) as f:
        counter = 0
        if pattern.match(file):
            print("Parsing {}.".format(file))
            totalsize = os.path.getsize(file)
            sumsize = 0
            for line in f:
                sumsize += sys.getsizeof(line)
                progress(sumsize, totalsize)
                if counter % 4 == 1:
                    splitline = line.strip().replace("\n", "")
                    splitline = splitline.split("N")
                    splitline = list(filter(None, splitline))
                    if splitline != []:
                        [sequences.append(k) for k in splitline]
                counter += 1
            return headers, sequences
        else:
            print("parsing {}".format(file))
            totalsize = os.path.getsize(file)
            sumsize = 0
            for line in f:
                sumsize += sys.getsizeof(line)
                progress(sumsize, totalsize)
                if line[0] == '>':
                    headers.append(line.strip())
                elif len(sequences) == len(headers):
                    splitline = line.replace('\n', '').split("N")   #remove \n and N from sequences
                    splitline = list(filter(None, splitline))       #remove empty strings (happens if multiple N's in a row)
                    if splitline == []:     #remove empty lists (happens if 'N\n' is the string)
                        continue
                    sequences[len(headers) - 1] += splitline.pop(0).strip()     #add rest of line to previous sequence
                    split_count = 1     #track headers for multiple sequences generated from one input sequence because of N's
                    for item in splitline:
                        headers.append("split_" + str(split_count))
                        split_count += 1
                        sequences.append(item)
                else:
                    sequences.append('')
                    splitline = line.replace('\n', '').split("N")
                    splitline = list(filter(None, splitline))
                    sequences[len(headers) - 1] += splitline.pop(0).strip()
                    split_count = 1
                    for item in splitline:
                        headers.append("split_" + str(split_count))
                        split_count += 1
                        sequences.append(item)
            return headers, sequences

#produces kmers from a list of sequences
def make_kmer_dict(seq_list):
    length = len(seq_list)
    kmer_dict = {}
    for index, seq in enumerate(seq_list):
        for character in range(len(seq) - kmer_size):
            hasher = hashlib.blake2b(bytes(seq[character : character + kmer_size], 'utf8'), digest_size = 7)
            kmer = int(hasher.hexdigest(), 16)
            if kmer in kmer_dict:
                kmer_dict[kmer] += 1
            else:
                kmer_dict[kmer] = 1
        if index % 10000 == 0:
            progress(index, length, "making kmers {:.0f} s elapsed".format(time.time() - starttime))
    progress(length, length, "making kmers {:.0f} s elapsed".format(time.time() - starttime))
    print()
    return kmer_dict

ref_headers, ref_sequences = split_headers(reference)
ref_hash = make_kmer_dict(ref_sequences)

D = sparse.csc_matrix((3,0))
for file in os.listdir(real_samples):
    headers, sequences = split_headers(real_samples + '/' + file)
    hashes = make_kmer_dict(sequences)
    intersection = set(ref_hash.keys()).intersection(input_hash.keys())
    intersection = tuple(intersection)
    interdist = [abs(ref_hash[x] - input_hash[x]) for x in intersection]
    distancedict = {str(intersection[x]) : str(interdist[x]) for x in range(len(intersection))}
    D.resize
    
