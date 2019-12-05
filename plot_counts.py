import sys
import matplotlib.pyplot as plt
import pandas as pd

'''USAGE: python plot_counts.py outfilename.pdf {'Tn-seq' or 'Bar-Seq'} filename1 filename2...'''

#import
outfilename=sys.argv[1]
input_type=sys.argv[2]
countfiles=sys.argv[3:]
countfile_ids=[]
for countfile in countfiles:
    countfile_ids.append(countfile.split('.')[0])


def get_BarSeq_counts(countfile):
    counts = []
    with open(countfile) as f:
        lenF=0.0
        for line in f:
            lenF+=1
            commsplit = line[2:-2].split(',')
            for i in commsplit:
                try:
                    n=float(i.split(': ')[1])
                    counts.append(normalized_n)
                except:
                    None
    f.close()
    
    normalized_counts = []
    for i in counts:
        normalized_counts.append(n/lenF)
        
    return normalized_counts

def get_TnSeq_counts(countfile):
    counts = []
    with open(countfile) as f:
        lenF=0.0
        for line in f:
            lenF+=1
            try:
                n=float(line.split('\t')[5])
                counts.append(n)
            except:
                None
    f.close()

##    normalized_counts = []
##    for i in counts:
##        normalized_counts.append(n/lenF)
    #print(normalized_counts)
    #return normalized_counts
    return counts

def plot_counts(all_counts):
    fig=plt.figure()
    ax = fig.add_subplot(1,1,1)

    ax.hist(all_counts, label=countfiles)
    plt.title('skew')
#    plt.xlabel('normalized abundance')
    plt.xlabel('abundance')
    plt.ylabel('number of inserts')
    plt.legend()
    plt.savefig(outfilename, bbox_inches='tight', format='pdf',dpi=1000)

def plot_low_counts(all_counts):
    fig=plt.figure()
    ax = fig.add_subplot(1,1,1)

    ax.hist(all_counts, label=countfiles)
    plt.title('skew')
#    plt.xlabel('normalized abundance')
    plt.xlabel('abundance')
    plt.ylabel('number of inserts')
    plt.xlim(0,1000)
    plt.legend()
    plt.savefig('lowcounts'+outfilename, bbox_inches='tight', format='pdf',dpi=1000)

def plot_all_counts():
    all_counts=[]
    if input_type=='Tn-Seq':
        for countfile in countfiles:
            all_counts.append(get_TnSeq_counts(countfile))
    elif input_type=='Bar-Seq':
        for countfile in countfiles:
            all_counts.append(get_BarSeq_counts(countfile))
    else:
        print('specify Tn-Seq or Bar-Seq as sys.argv[1]')
    plot_counts(all_counts)
    plot_low_counts(all_counts)

plot_all_counts()
