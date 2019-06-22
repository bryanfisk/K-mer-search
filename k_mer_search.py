print("importing packages")
from scipy import sparse
import hashlib
import argparse
import time
import sys
import os
import random
random.seed(42)
start_time = time.time()
kmer_size = 15
file_dir = "D:/Desktop/metagenomes/test_folder"#finished"
reference = "D:/Desktop/metagenomes/clostridioides_difficile_genome.fna"

#get options from command line input
parser = argparse.ArgumentParser(description = 'Generate reverse complement for sequence input.')
parser.add_argument('--output', '-o', help = 'Set output file.')
parser.add_argument('--input', '-i', help = 'Set input file.')
parser.add_argument('--kmer', '-k', const = 10, type = int, help = 'Set k-mer size. (10)', nargs = '?')
args = parser.parse_args()

#progress bar
def progress(count, total, status = ''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))
    percents = int(100.0 * count / float(total))
    bar = '■' * filled_len + '-' * (bar_len - filled_len)
    sys.stdout.write('|{}| {}% {}\r'.format(bar, percents, status))
    sys.stdout.flush()

def format_time(start, end):
    elapsed = end - start
    hours = int(elapsed // 3600)
    minutes = int(elapsed - hours * 3600) // 60
    seconds = int(elapsed - hours * 3600 - minutes * 60)
    return '{:02d}:{:02d}:{:02d}'.format(hours, minutes, seconds)

#split sequences from headers in fastq file
def split_headers(file):
    sequences = []
    with open(file) as f:
        if '.fastq' in file:
            counter = 0
            print("Parsing {}.".format(file))
            totalsize = os.path.getsize(file)
            sumsize = 0
            for line in f:
                sumsize += len(line)
                if sumsize % 30000 == 0 or sumsize == totalsize:
                    progress(sumsize, totalsize, 'parsing')
                if counter % 4 == 1:
                    splitline = line.strip().replace("\n", "")
                    splitline = splitline.split("N")
                    splitline = list(filter(None, splitline))
                    if splitline != []:
                        [sequences.append(k) for k in splitline]
                counter += 1
            return sequences
        else:
            print("parsing {}.".format(file))
            totalsize = os.path.getsize(file)
            sumsize = 0
            cat_line = ''
            for line in f:
                sumsize += len(line)
                if sumsize % 10000 == 0 or sumsize == totalsize:
                    progress(sumsize, totalsize, 'parsing')
                if line[0] == '>':
                    sequences, cat_line = clean_sequence(sequences, cat_line)
                else:
                    cat_line = cat_line + line
            sequences, cat_line = clean_sequence(sequences, cat_line)
            return sequences

def clean_sequence(sequence_list, cat):
    if cat != '':
        cat = cat.replace('\n', '').split('N')
        cat = list(filter(None, cat))
        [sequence_list.append(k) for k in cat]
    cat = ''
    return sequence_list, cat

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
            time_elapsed = time.time() - start_time
            progress(index, length, "making kmers {} elapsed".format(format_time(start_time, time.time())))
    progress(length, length, "making kmers {} elapsed".format(format_time(start_time, time.time())))
    print()
    return kmer_dict

#take directory and reference input and return matrix of distance values of kmers in common
def dir_to_matrix(dir, ref):
    #set start time, parse reference, hash reference
    ref_seq = split_headers(ref)
    ref_hash = make_kmer_dict(ref_seq)
    X = sparse.csc_matrix((0, 0))
    y = sparse.csc_matrix((0, 0))
    header_dict = {}
    new_item_count = 0
    for index, file in enumerate(os.listdir(dir)):
        #columns = X.shape[1]
        
        #make label matrix
        if 'pos_' in file:
            y = sparse.hstack((y, [1]))
        else:
            y = sparse.hstack((y, [0]))
        
        #parse input, hash input, get input-reference common kmer distances
        input_seq = split_headers(file_dir + '/' + file)
        input_hash = make_kmer_dict(input_seq)
        intersection = tuple(set(ref_hash.keys()).intersection(input_hash.keys()))
        for item in intersection:
            if item not in header_dict:
                header_dict[item] = new_item_count
                new_item_count += 1
        interdist = [abs(ref_hash[x] - input_hash[x]) for x in intersection]
        rows = [0 for _ in range(len(interdist))]
        intersection = [header_dict[x] for x in intersection]
        temp_matrix = sparse.csc_matrix((interdist, (rows, intersection)))
        
        #pad either matrix depending on which is shorter
        if temp_matrix.shape[1] > X.shape[1]:#columns:
            X.resize((X.shape[0], temp_matrix.shape[1]))
        elif temp_matrix.shape[1] < X.shape[1]:#columns:
            temp_matrix.resize((1, X.shape[1]))

        #append temp_matrix to X
        X = sparse.vstack((X, temp_matrix))
    return X, y

#reduce dimensionality by removing columns with the same values
def remove_redundant_columns(D):
    D = D.tocsc()
    start_length = D.shape[1]
    print('starting length:', start_length)
    keepset = [x for x in range(D.shape[1])]
    for column in range(D.shape[1]):
        if column % 10000 == 0:
            progress(column, D.shape[1], "reducing dims, {} elapsed".format(format_time(start_time, time.time())))
        if max(D[:,column].toarray()) == min(D[:,column].toarray()):
            keepset.remove(column)
    print(D.shape)
    D = D[:, keepset]
    print('\n' ,start_length, "columns reduced to", D.shape)
    return D

X, y = dir_to_matrix(file_dir, reference)
X = remove_redundant_columns(X)

#sparse.save_npz('X.npz', X)
#sparse.save_npz('y.npz', y)