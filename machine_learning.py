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
real_samples = "C:/Users/bryan/Desktop/metagenomes/reallifesamples"
starttime = time.time()

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
    #os.system('cls')
    pattern = re.compile('.*\.fastq$')
    headers = []
    sequences = []
    with open(file) as f:
        counter = 0
        if pattern.match(file):
            print("Parsing {}.".format(file))
            for line in f:
                if counter % 4 == 1:
                    splitline = line.strip().replace("\n", "")
                    splitline = splitline.split("N")
                    splitline = list(filter(None, splitline))
                    if splitline != []:
                        [sequences.append(k) for k in splitline]
                    #sequences.append(line.strip().replace("N", ""))
                counter += 1
            return headers, make_kmer_dict(sequences, file)
        else:
            print("parsing {}".format(file))
            for line in f:
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
            return headers, make_kmer_dict(sequences, file)

#produces kmers from a list of sequences
def make_kmer_dict(seq_list, file):
    length = len(seq_list)
    #if length > 50:
    #    twopercent = length // 50
    #else:
    #    twopercent = 1
    kmer_dict = {}
    for index, seq in enumerate(seq_list):
        #totalpercent = math.ceil((index + 1) / length * 100)
        #if index % twopercent == 0:
            #os.system('cls')
            #blocks = math.ceil((index + 1) / length * 50)
            #noblocks = 50 - blocks
            #print("Making k-mers from {}.".format(file))
            #print("█" * blocks, "-" * noblocks, ' {:d}%'.format(totalpercent), sep = "")
        for character in range(len(seq) - kmer_size):
            hasher = hashlib.blake2b(bytes(seq[character : character + kmer_size], 'utf8'), digest_size = 7)
            kmer = int(hasher.hexdigest(), 16)
            if kmer in kmer_dict:
                kmer_dict[kmer] += 1
            else:
                kmer_dict[kmer] = 1
        if index % 10000 == 0:
            progress(index, length, "making kmers {:.0f} seconds elapsed".format(time.time() - starttime))
    progress(length, length, "making kmers {:.0f} seconds elapsed".format(time.time() - starttime))
    print()
    return kmer_dict

for file in os.listdir(real_samples):
    headers, hashes = split_headers(real_samples + '/' + file)
