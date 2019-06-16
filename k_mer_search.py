print("importing packages")
#import matplotlib.pyplot as plt
import hashlib
import argparse
import time
import math
import sys
import os
import re

kmer_size = 15
file_dir = "C:/Users/bryan/Desktop/metagenomes/finished/"
#input = "ecoli.fastq"
#reference = "e.coli_genome.fna"
reference = "C:/Users/bryan/Desktop/metagenomes/clostridioides_difficile_genome.fna"

#get files from root directory
for _, _, files in os.walk(file_dir):
    pass

import random
import string
random.seed(42)	

def numsplit(string, n):
	if len(string) < n: return string
	numsplits = len(string) // n
	return [string[x : x + n] for x in range(0, numsplits * n, n)]
		
def pad(string, n):
	if len(string) < n:
		numzeroes = n - len(string)
		string = string + [0 for _ in range(numzeroes)]
	elif len(string) > n:
		totallength = (len(string) // n) * n + n
		numzeroes = totallength - len(string)
		string = string + [0 for _ in range(numzeroes)]
	return string

#progress bar
def progress(count, total, status = ''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))
    percents = round( 100.0 * count / float(total), 1)
    bar = '■' * filled_len + '-' * (bar_len - filled_len)
    sys.stdout.write('|{:}| {:}{:}  {:}\r'.format(bar, percents, '%', status))
    sys.stdout.flush()

#get options from command line input
parser = argparse.ArgumentParser(description = 'Generate reverse complement for sequence input.')
parser.add_argument('--output', '-o', help = 'Set output file.')
parser.add_argument('--input', '-i', help = 'Set input file.')
parser.add_argument('--kmer', '-k', const = 10, type = int, help = 'Set k-mer size. (10)', nargs = '?')
args = parser.parse_args()

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

#Takes list of tuples and combines tuples that overlap one another. Leaves non-overlapping tuples alone.
def combine(t):
    t.sort()
    output = []
    for tup in t:
        if len(output) == 0:
            output.append(tup)
        else:
            flag = True
            for i, o in enumerate(output):
                if tup[0] < o[0] and tup[1] >= o[0]:
                    output[i] = (tup[0], o[1])
                    flag = False
                if tup[0] <= o[1] and tup[1] > o[1]:
                    output[i] = (o[0], tup[1])
                    flag = False
            if (tup[0] > o[1] or tup[1] < o[0]) and flag == True:
                output.append(tup)
    return output

def compile_uniques(dictionary, intersection):
    common = [dictionary[x] for x in intersection]
    common.sort()
    common_sum = sum(common)
    common = [x / common_sum for x in common]
    return common

def find_indices(d, s):
    temp = []
    for x in range(len(s)):
        if s[x] in d:
            temp.append((x, d[s[x]]))
    return temp

starttime = time.time()

ref_headers, ref_hash = split_headers(reference)

#Make dictionary with each hashsed kmer as a key 
#and the number of times that kmer is encountered as a value     

ref_unique_kmers = [k for k, v in ref_hash.items() if v == 1]

fileout = open('output.txt', 'w+')
masterset = set()
feature_dict_raw = {}
for index, file in enumerate(files):
    input_headers, input_hash = split_headers(file_dir + file)
    input_unique_kmers = [k for k, v in input_hash.items() if v == 1]

    intersection = set(ref_hash.keys()).intersection(input_hash.keys())
    intersection = tuple(intersection)
    interdist = [abs(ref_hash[x] - input_hash[x]) for x in intersection]

    #ref_common = compile_uniques(ref_hash, intersection)
    #input_common = compile_uniques(input_hash, intersection)

    #diff = sum([abs(input_common[x] - ref_common[x]) for x in range(len(input_common))])
    #score = round(-math.log(diff/2), 2)
    #jaccard_dist = len(input_common) / (len(ref_unique_kmers) + len(input_unique_kmers) + len(input_common))

    distancedict = {str(intersection[x]) : str(interdist[x]) for x in range(len(intersection))}

    keys = list(distancedict.keys())
    values = list(distancedict.values())

    fileout.write(file + '\n')
    fileout.write(','.join(keys) + '\n')
    fileout.write(','.join(values) + '\n')
    #feature_dict_raw[file] = distancedict

'''
for i in intersection:
    if i not in masterset:
        masterset.add(i)

masterset = tuple(masterset)
feature_dict_final = {}
count = 1
for key, value in list(feature_dict_raw.items()):
    feature_dict_final[key] = ['0' for _ in masterset]
    index_value_pairs = find_indices(value, masterset)
    for index, value in index_value_pairs:
        if count % 1000 == 0:
            print(count)
        feature_dict_final[key][index] = value
        count += 1

with open("output.csv", 'w') as output:
    output.write('label,' + ','.join(list(map(str, masterset))) + '\n')
    for label, feature in list(feature_dict_final.items()):
        temp = str(feature)
        output.write(label + ',' + ','.join(feature))



    #plot stuff
    fig = plt.figure()
    ax1 = fig.add_axes((0.15, 0.25, 0.8, 0.7))
    ax1.plot(ref_common)
    ax1.plot(input_common)
    ax1.set_title('Frequencies of Shared k-mers')
    ax1.set_ylabel('K-mer frequency', rotation = 90)
    ax1.set_xlabel('K-mer index')
    ax1.legend(loc = 'upper left', labels = ('Reference', 'Input'))
    fig.text(0.1, 
             0.03, 
             'Different reference k-mers: {}     Unique reference k-mers: {}\n'.format(len(ref_hash), len(ref_unique_kmers)), ha = 'left')
    fig.text(0.1,
             0.065, 
             'Different input k-mers: {}            Unique input k-mers: {}\n'.format(len(input_hash), len(input_unique_kmers)), ha = 'left')
    fig.text(0.1,
             0.1,
             'Number of shared k-mers: {}       Jaccard Distance: {:.2f}\n'.format(len(intersection), jaccard_dist), ha = 'left')
    plt.show()
    '''