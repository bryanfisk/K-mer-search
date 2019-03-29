import argparse
import re

k_mer_size = 15

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
        print("Parsing fastq input.")
        if pattern.match(file):
            for line in f:
                if counter % 4 == 1:
                    sequences.append(line.strip())
                counter += 1
            return headers, make_kmers(sequences)
        else:
            print("Parsing fna reference.")
            for line in f:
                if line[0] == '>':
                    headers.append(line.strip())
                elif len(sequences) == len(headers):
                    sequences[len(headers) - 1] += line.strip()
                else:
                    sequences.append('')
                    sequences[len(headers) - 1] += line.strip()
            return headers, make_kmers(sequences)

#produces kmers from a list of sequences
def make_kmers(seq_list):
    print("Making kmers.")
    length = len(seq_list)
    fivepercent = length * 5 // 100
    kmer_list = []
    count = 0
    for seq in seq_list:
        if count % fivepercent == 0 and count > fivepercent - 1:
            print(round(count * 100 / length), "% done.")
        kmer_sublist = []
        for character in range(len(seq) - k_mer_size):
            kmer_sublist.append(seq[character : character + k_mer_size])
        kmer_list.append(kmer_sublist)
        count += 1
    print("Hashing kmers.")
    kmer_hash = [[hash(kmer) for kmer in sub_list] for sub_list in kmer_list]
    return kmer_hash

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

hash_dict = {hash:"" for sequence in input_hash for hash in sequence}

nc = 0
match_count = 0
print("Comparing kmers.")
for contig in ref_hash:
    nc += 1 
    print("contig", nc, "of", len(ref_hash))
    for hash in contig:
        if hash in hash_dict:
            match_count += 1

print(match_count)
total_length = sum([len(k) for k in ref_hash])
print(total_length)
percent_match = match_count / total_length

print('Percent reference genome kmers matching input kmers: ', percent_match )