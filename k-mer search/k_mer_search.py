k_mer_size = 5

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

test = [(1,3), (5,7), (2,6), (11,14)]
print(combine(test))

#Replacing below with k-1-mer hashing and lookup in pathogen genome library
'''
headers = []
sequences = []
with open("metagenome_sample.txt") as file:
    for line in file:
        if line[0] == '>':
            header.append([line, [None]])
        else:
            sequences.append(line)
    file.close()

k_mer_list = []
with open("pathogen_sequences.txt") as file:
    for line in file:
        for c in len(line - k_mer_size):
            k_mer_list.append(line[c : c + k_mer_size])
    file.close()

for k_mer in k_mer_list:
    for index, line in enumerate(sequences):
        found = line.find(k_mer)
        if found and header[index][1][0] == None:
            header[index][1][0] = (found, found + )
'''