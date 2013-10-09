#!/usr/bin/env python
"""
Write a JTK that does the following:
1) Takes in a list of ZTs and gene values
2) Allows a choice of waveform
3) Allows a choice of phase and period
4) Takes the number of available points for a given gene and calculates the null distribution for that set of timepoints
5) Calculates the Kendall's tau between the time points and the Null distribution
"""

from scipy.stats import kendalltau
import numpy as np
import sys

def main():
    fn=sys.argv[1]
    read_in(fn)


def read_in(fn):
    with open(fn,'r') as f:
        master_match=[]
        for line in f:
            #print line
            words = line.strip().split()
            words=[word.strip() for word in words]
            #print words
            match_searched = []
            if words[0] =="#":
                print words
                for i in xrange(1,len(words)):
                    match_searched.append(i)
                    match = [i]
                    for j in xrange(i+1,len(words)):
                        if j not in match_searched:
                            if words[i].split('_')[0] == words[j].split('_')[0]:
                                match.append(j)
                                match_searched.append(j)
                                    
                    if len(match) > 1:
                        master_match.append(match)
                print master_match
            else:
                new = [words[0]]
                #print set([m for match in master_match for m in match])
                to_collapse = set([m for match in master_match for m in match])
                to_check = set(range(1,len(words)))
                for i in range(1,len(words)):
                    if i in to_check:
                        #print i,"is i here"
                        #print to_collapse,"are indicies left"
                        if i in to_collapse:
                            for match in master_match:
                                if i in match:
                                    sum = 0
                                    for m in match:
                                        sum += float(words[m])
                                    sum = sum / len(match)
                                    new.append(sum)
                                    to_collapse = to_collapse.difference(match)
                                    to_check = to_check.difference(match)
                        else:
                            new.append(float(words[i]))
                print "new is",new


def combine_replicates():
    pass
"""
If header has biological replicates, combine them in all lists, return list of lists
"""    



def organize_data():
    pass
"""
Organize list of lists from such that genes with similar time-series holes match (for null distribution calc)
"""






if __name__=="__main__":
    main()
