import matplotlib.pyplot as plt
import argparse
import math
import os
import re

kmer_size = 4

#get options from command line input
parser = argparse.ArgumentParser(description = 'Generate reverse complement for sequence input.')
parser.add_argument('--output', '-o', help = 'Set output file.')
parser.add_argument('--input', '-i', help = 'Set input file.')
parser.add_argument('--kmer', '-k', const = 10, type = int, help = 'Set k-mer size. (10)', nargs = '?')
args = parser.parse_args()

#split headers from sequences from fastq file
def split_headers(file):
    os.system('cls')
    pattern = re.compile('.*\.fastq$')
    headers = []
    sequences = []
    with open(file) as f:
        counter = 0
        if pattern.match(file):
            print("Parsing {}.".format(file))
            for line in f:
                if counter % 4 == 1:
                    sequences.append(line.strip().replace("N", ""))
                counter += 1
            return headers, make_kmer_dict(sequences, file)
        else:
            print("Parsing {}.".format(file))
            for line in f:
                if line[0] == '>':
                    headers.append(line.strip())
                elif len(sequences) == len(headers):
                    sequences[len(headers) - 1] += line.strip()
                else:
                    sequences.append('')
                    sequences[len(headers) - 1] += line.strip()
            return headers, make_kmer_dict(sequences, file)

#produces kmers from a list of sequences
def make_kmer_dict(seq_list, file):
    length = len(seq_list)
    if length > 1:
        twopercent = length // 50
    else:
        twopercent = 1
    kmer_dict = {}
    for index, seq in enumerate(seq_list):
        totalpercent = math.ceil((index + 1) / length * 100)
        if index % twopercent == 0:
            os.system('cls')
            blocks = math.ceil((index + 1) / length * 50)
            noblocks = 50 - blocks
            print("Making k-mers from {}.".format(file))
            print("█" * blocks, "-" * noblocks, ' {:d}%'.format(totalpercent), sep = "")
        for character in range(len(seq) - kmer_size):
            kmer = hash(seq[character : character + kmer_size])
            if kmer in kmer_dict:
                kmer_dict[kmer] += 1
            else:
                kmer_dict[kmer] = 1
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

input = "filtered_SRR061908.fastq"
reference = "betapapillomavirus.fna"

input_headers, input_hash = split_headers(input)#args.input)
ref_headers, ref_hash = split_headers(reference)#reference.txt)

#Make dictionary with each hashsed kmer as a key 
#and the number of times that kmer is encountered as a value     

ref_unique_kmers = [k for k, v in ref_hash.items() if v == 1]
input_unique_kmers = [k for k, v in input_hash.items() if v == 1]

print('Different reference k-mers: {}\nReference unique k-mers: {}\nDifferent input k-mers: {}\nInput unique k-mers: {}'.format(len(ref_hash), len(ref_unique_kmers), len(input_hash), len(input_unique_kmers)))

print('Intersection')
intersection = set(ref_hash.keys()).intersection(input_hash.keys())
print('Number of shared k-mers:', len(intersection))

print('Compiling and normalizing unique reference k-mers counts.')
ref_common = compile_uniques(ref_hash, intersection)

print('Compiling and normalizing unique input k-mers counts.')
input_common = compile_uniques(input_hash, intersection)

print('Calculating similarity score.')
diff = sum([abs(input_common[x] - ref_common[x]) for x in range(len(input_common))])
score = round(-math.log(diff/2), 2)

print('Similarity:', score)

plt.plot(ref_common)
plt.plot(input_common)
plt.title('Frequencies of Shared k-mers')
plt.legend(loc = 'upper left', labels = ('Reference', 'Input'))
plt.show()