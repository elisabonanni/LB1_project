''' Extract sequences from a fasta file taking as input the ids of the sequences.'''

#!/usr/bin/env python3
import sys

def get_seqs(fileseq,pos=1):
    '''Read a fasta file and create a dictionary with keys=ids and values=corresponding sequences'''
    dseq={} #initialise the dictionary
    with open(fileseq, "r") as f: #read the file
        for line in f:
            if line[0]==">":
                k = line[1:].rstrip().split("|")[pos] #starts from line[1:] to exclude the > sign
                dseq[k]='' #initialise the dictionary sequence
                continue #go to the next line
            dseq[k] = dseq[k] + line.rstrip()
    return dseq

if __name__=="__main__":
    fileseq = sys.argv[1]
    fileids = sys.argv[2]
    pos=1 #standardly the id is in the first positition
    if len(sys.argv)>3:
        pos=int(sys.argv[3])-1 #if exists an extra input, it corresponds to the position of the ids
    dseq = get_seqs(fileseq,pos)
    ids=open(fileids, "r").read().rstrip().split('\n') #create a list of the identifiers from the file
    for i in ids:
        if dseq.get(i,0) != 0: #if the id is not present
            print(">"+i+"\n"+dseq[i])
        else:
            print('WARNING: Sequence '+i+' not found', file=sys.stderr)
