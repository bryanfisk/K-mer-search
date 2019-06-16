print("importing packages")
import sys
import os
import re
from scipy import sparse
import numpy as np
import time

masterdict = {}
dictcounter = 0
y_pattern = re.compile('^pos_')
y = []
count = 0
line_count = 0
D = sparse.csc_matrix((952,0))
elapsedtime = 0

def progress(count, total, status = ''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))
    percents = round( 100.0 * count / float(total), 1)
    bar = '■' * filled_len + '-' * (bar_len - filled_len)
    sys.stdout.write('|{:}| {:}{:}  {:}\r'.format(bar, percents, '%', status))
    sys.stdout.flush()


def neg_to_pos_hash(hash):
    if hash < 0:
        return hash % (sys.maxsize)
    else:
        return hash

#read input
totalsize = os.path.getsize('output.txt')
bytecount = 0
print("reading file")
with open("output.txt", 'r') as file:
    starttime = time.time()
    for line in file:
        bytecount += sys.getsizeof(line)
        #first line in each triplet of lines, get y values from regex matching 'pos' in sample name
        if line_count % 3 == 0:
            label = line.strip()
            if y_pattern.search(label):
                y.append(1)
            else:
                y.append(0)
        #second line in each triplet of lines, get hash keys
        elif line_count % 3 == 1:
            splitkeys = line.strip()
            splitkeys = splitkeys.split(',')
            splitkeys = map(int, splitkeys)
            splitkeys = list(splitkeys)
            splitkeys = [neg_to_pos_hash(k) for k in splitkeys]
        #third line in each triplet of lines, get matrix values
        else:
            splitvalues = line.strip()
            splitvalues = splitvalues.split(',')
            splitvalues = map(int, splitvalues)
            splitvalues = list(splitvalues)
            
            #add data to matrix
            for key in splitkeys:
                if key not in masterdict.keys():
                    masterdict[key] = dictcounter
                    dictcounter += 1
            splitkeys = [masterdict[k] for k in splitkeys]
            rownums = [count for _ in range(len(splitkeys))]
            tempmatrix = sparse.csc_matrix((splitvalues, (rownums, splitkeys)), shape = (952, dictcounter))
            tempmatrix.resize(952, dictcounter)
            D.resize(952, dictcounter)
            D = D + tempmatrix
            elapsedtime = time.time() - starttime
            progress(bytecount, totalsize, "sample {:d}, {:d} columns, last sample: {:.2f} sec".format(count + 1, D.shape[1], elapsedtime))
            count += 1
        line_count += 1

y = np.asarray(y)
np.savez("y_values.npz", y)

#feature engineering

print()
startlength = D.shape[1]
starttime = time.time()
keepset = [x for x in range(D.shape[1])]
for column in range(D.shape[1] - 1):
    if column % 10000 == 0:
        progress(column, D.shape[1], "removing redundant columns, elapsed time: {:.2f}".format(time.time() - starttime))
    if max(D[:,column].toarray()) == min(D[:,column].toarray()):
        keepset.remove(column)
D = D[:, keepset]
print('\n' ,startlength, "columns reduced to", D.shape[1])

sparse.save_npz("extracted_features.npz", D)
