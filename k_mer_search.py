import argparse
import math
import os
import re

kmer_size = 15

#get options from command line input
parser = argparse.ArgumentParser(description = 'Generate reverse complement for sequence input.')
parser.add_argument('--output', '-o', help = 'Set output file.')
parser.add_argument('--input', '-i', help = 'Set input file.')
parser.add_argument('--kmer', '-k', const = 10, type = int, help = 'Set k-mer size. (10)', nargs = '?')
args = parser.parse_args()

#split headers from sequences from fastq file
def split_headers(file):
    pattern = re.compile('.*\.fastq$')
    headers = []
    sequences = []
    with open(file) as f:
        counter = 0
        print("Parsing fastq.")
        if pattern.match(file):
            for line in f:
                if counter % 4 == 1:
                    sequences.append(line.strip())
                counter += 1
            return headers, make_kmer_dict(sequences, file)
        else:
            print("Parsing fna.")
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
                pass
            else:
                kmer_dict[kmer] = ""
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

input = "ecoli.fastq"
reference = "E. coli genome.fna"

input_headers, input_hash = split_headers(input)#args.input)
ref_headers, ref_hash = split_headers(reference)#reference.txt)

#hash_dict = {hash:"" for sequence in input_hash for hash in sequence}

match_count = 0
print("Comparing kmers.")
for ref_kmer in ref_hash:
    if ref_kmer in input_hash:
        match_count += 1

percent_match = match_count / len(ref_hash) * 100
prob_one_match = len(ref_hash) / 4 ** kmer_size
average_match = prob_one_match * len(input_hash)
average_prop_match = average_match / match_count
probability = ((4 ** kmer_size) - len(ref_hash)) / (4 ** kmer_size)

print('Percent reference genome kmers matching input kmers: {:.2f}%'.format(percent_match))
print('Average expected number of random matches: {:f}\nproportion of all matches: {:.2f}'.format(average_match, average_prop_match))

