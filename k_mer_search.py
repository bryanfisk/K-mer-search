print("importing packages")
import argparse
from functions import dir_to_matrix
from scipy import sparse
import csv

file_dir = "D:/Desktop/metagenomes/finished"
reference = "D:/Desktop/metagenomes/clostridioides_difficile_genome.fna"
sample_files = "D:/Desktop/metagenomes/test_folder"

#get options from command line input
parser = argparse.ArgumentParser(description = 'Generate reverse complement for sequence input.')
parser.add_argument('--output', '-o', help = 'Set output file.')
parser.add_argument('--input', '-i', help = 'Set input file.')
parser.add_argument('--kmer', '-k', const = 10, type = int, help = 'Set k-mer size. (10)', nargs = '?')
args = parser.parse_args()

X, y, index_dictionary = dir_to_matrix(sample_files, reference)

sparse.save_npz('X.npz', X)
sparse.save_npz('y.npz', y)
with open("index_dictionary.tsv", "w") as file:
    for key, value in index_dictionary.items():
        file.write('{}\t{}\n'.format(str(key), str(value)))