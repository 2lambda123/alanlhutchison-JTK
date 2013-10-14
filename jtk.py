#!/usr/bin/env python
"""
You should probably be calling this from wrap_jtk.sh

Write a JTK that does the following:
1) Takes in a list of ZTs and gene values
2) Allows a choice of waveform
3) Allows a choice of phase and period
4) Takes the number of available points for a given gene and calculates the null distribution for that set of timepoints
5) Calculates the Kendall's tau between the time points and the Null distribution
"""
VERSION="0.0"
from scipy.stats import kendalltau
from operator import itemgetter
import numpy as np
import sys
import argparse

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.stats import norm


def main(args):
    fn = args.filename
    waveform = args.function
    period = args.period
    phase = args.phase
    eff_amp = 1
    updated = read_in(fn)
    header,series =organize_data(updated)
    reference = generate_base_reference(header,waveform,phase,period)
    #print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format("#ID","waveform","period","phase","eff_amp","tau","p")
    RealKen = KendallTauP()        
    for serie in series:
        geneID,tau,p = generate_mod_series(reference,serie,RealKen)
        print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(geneID,waveform,period,phase,eff_amp,tau,p)


def read_in(fn):
    """Reads in a file in correct '#\tZTX_X\tZTX_X\tZTX_X\n geneID\tvalue\tvalue'
       Returns a list of lists of lines with replicates combined
       Correctly deals with the NA case:
           If NA is part of a replicate, ignore from mean
           If NA is alone or all replicates, propagate through to output"""
    updated = []
    with open(fn,'r') as f:
        master_match=[]
        for line in f:
            words = line.strip().split()
            words=[word.strip() for word in words]
            if words[0] =="#":
                match_searched = []
                #print words
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
                #print master_match

            new = [words[0]]
            to_collapse = set([m for match in master_match for m in match])
            to_check = set(range(1,len(words)))
            for i in range(1,len(words)):
                if i in to_check:
                    if i in to_collapse:
                        for match in master_match:
                            if i in match:
                                sum = 0
                                if "ZT" in line:
                                    sum = words[i].split('_')[0]
                                else:
                                    NAflag = 0
                                    for m in match:
                                        if words[m].strip()!="NA":
                                            sum += float(words[m])
                                        else:
                                            NAflag += 1
                                    if NAflag != len(match):
                                        sum = sum / float(len(match)-NAflag)
                                    elif NAflag==len(match): 
                                        sum = "NA"
                                new.append(sum)
                                to_collapse = to_collapse.difference(match)
                                to_check = to_check.difference(match)
                    else:
                        if "ZT" in line:
                            new.append(words[i].split('_')[0])
                        else:
                            if words[i]!="NA":
                                new.append(float(words[i]))
                            else:
                                new.append(words[i])
            #print "new is",new
            updated.append(new)
#    for update in updated:
#        print update
    return updated




def organize_data(updated):
    """
    Organize list of lists from such that genes with similar time-series holes match (for null distribution calc)
    Return a header ['#','ZTX','ZTY'...] and a list of lists [ lists with similar holes (identical null distribution) , [],[],[]] 
    """
    header = updated[0]
    L = updated[1:]

    for i in xrange(1,len(header)):
        L=sorted(L, key=itemgetter(i))

#    print "Header is"
#    print header
#    for line in L:
#        print line
    return header,L

def generate_base_reference(header,waveform="cosine",phase=0,period=24):
    """
    This will generate a waveform with a given phase and period based on the header, 
    """
    tpoints = []
    ZTs = header[1:]
    coef = 2.0 * np.pi / float(period)
    for ZT in ZTs:
        z = ZT[2:]
        tpoints.append( (float(z)+float(phase) ) * coef)
    #print tpoints
    #print [tpoint/np.pi/2.0 for tpoint in tpoints]

    if waveform == "cosine":
        reference=np.cos(tpoints)
    elif waveform == "impulse":
        reference=np.cos(tpoints)        
    elif waveform == "rampup":
        reference=np.cos(tpoints)
    elif waveform == "rampdown":
        reference=np.cos(tpoints)
    elif waveform == "step":
        reference=np.cos(tpoints)
    

    return reference


def generate_mod_series(reference,series,RealKen):
    """
    Takes the series from generate_base_null, takes the list from data, and makes a null
    for each gene in data or uses the one previously calculated.
    Then it runs Kendall's Tau on the exp. series against the null
    """
    geneID = series[0]
    values = series[1:]
    binary = [1 if value!="NA" else np.nan for value in values]
    temp = reference*binary
    mod_reference = [value for value in temp if not np.isnan(value)]
    mod_values = [value for value in values if value!='NA']
#    print reference
#    print temp
#    print mod_reference
#    print mod_values

    if len(mod_values) < 3:
        tau,p = np.nan,np.nan
    elif mod_values.count(np.nan) == len(mod_values):
        tau,p = np.nan,np.nan
    elif mod_values.count(0) == len(mod_values):
        tau,p = np.nan,np.nan
    elif sum(mod_values)<0.00001:
        tau,p = np.nan,np.nan        
    else:
        tau,p=kendalltau(mod_values,mod_reference)
        if not np.isnan(tau):
            pk = RealKen.pval(tau,len(mod_values))
            if pk!=None:
                p=pk
    
    #print tau,p
    return geneID,tau,p




def __create_parser__():
    p = argparse.ArgumentParser(
        description="python script runner for JTK_CYCLE statistical test",
        epilog="...",
        version=VERSION
        )

                   
    p.add_argument("-t", "--test",
                   action='store_true',
                   default=False,
                   help="run the Python unittest testing suite")

    p.add_argument("-f", "--file",
                   dest="filename",
                   action='store',
                   metavar="FILENM",
                   type=str,
                   help="give a filename else this thang won't run")

    analysis = p.add_argument_group(title="JTK_CYCLE analysis options")
    analysis.add_argument("--function",
                          dest="function",
                          type=str,
                          metavar="$FUNC_STR",
                          action='store',
                          default="cosine",
                          choices=["cosine","rampup","rampdown","step","impulse"],
                          help="cosine (dflt), rampup, rampdown, impulse, step")
    analysis.add_argument("-w", "--width",
                          dest="width",
                          type=float,
                          metavar="W",
                          action='store',
                          default=0.75,
                          help="shape parameter for alt. waveforms \in [0,1]")
    analysis.add_argument("-ph", "--phase",
                          dest="phase",
                          metavar="P",
                          type=float,
                          default=0.0,
                          help="set phase of reference waveform (dflt: 0.0)")
    analysis.add_argument("-p","--period",
                         dest="period",
                         metavar=float,
                         type=float,
                         action='store',
                         help="set period to be searched")

    
    distribution = analysis.add_mutually_exclusive_group(required=False)
    distribution.add_argument("-e", "--exact",
                              dest="harding",
                              action='store_true',
                              default=False,
                              help="use Harding's exact null distribution (dflt)")
    distribution.add_argument("-n", "--normal",
                              dest="normal",
                              action='store_true',
                              default=False,
                              help="use normal approximation to null distribution")
    
    
    return p


# instantiate class to precalculate distribution
# usage: 
#   K = KendallTauP()
#   pval = K.pval(tau,n,two_tailed=True)
class KendallTauP:
    def __init__(self,N=25):        
        # largest number of samples to precompute
        self.N = N
        Nint = self.N*(self.N-1)/2

        # first allocate freq slots for largest sample array
        # as we fill this in we'll save the results for smaller samples

        # total possible number of inversions is Nint + 1
        freqN = np.zeros(Nint + 1)
        freqN[0] = 1.0

        # save results at each step in freqs array
        self.freqs = [np.array([1.0])]
        for i in xrange(1,self.N):
            last = np.copy(freqN)
            for j in xrange(Nint+1):
                # update each entry by summing over i entries to the left
                freqN[j] += sum(last[max(0,j-i):j])
            # copy current state into freqs array
            # the kth entry of freqs should have 1+k*(k-1)/2 entries
            self.freqs.append(np.copy(freqN[0:(1+(i+1)*i/2)]))
            
        # turn freqs into cdfs
        # distributions still with respect to number of inversions
        self.cdfs = []
        for i in xrange(self.N):
            self.cdfs.append(np.copy(self.freqs[i]))
            # turn into cumulative frequencies
            for j in xrange(1,len(self.freqs[i])):
                self.cdfs[i][j] += self.cdfs[i][j-1]
            # convert freqs to probs
            self.cdfs[i] = self.cdfs[i]/sum(self.freqs[i])
            
    # plot exact distribution compared to normal approx
    def plot(self,nlist):
        colors = cm.Set1(np.linspace(0,1,len(nlist)))

        # for plotting gaussian
        x = np.linspace(-1.2,1.2,300)
        # plot pdfs
        plt.figure()
        for i in xrange(len(nlist)):
            ntot = len(self.freqs[nlist[i]-1])-1
            tauvals = (ntot - 2.0*np.arange(len(self.freqs[nlist[i]-1])))/ntot
            probs = ((ntot+1.0)/2.0)*self.freqs[nlist[i]-1]/sum(self.freqs[nlist[i]-1])
            plt.scatter(tauvals,probs,color=colors[i])
            # now plot gaussian comparison
            var = 2.0*(2.0*nlist[i]+5.0)/(nlist[i]*(nlist[i]-1)*9.0)
            plt.plot(x,norm.pdf(x,0.0,np.sqrt(var)),color=colors[i])
        plt.legend(nlist,loc='best')
        # plt.savefig('pdfs.png')
        plt.show()

        # now plot cdfs
        plt.figure()
        for i in xrange(len(nlist)):
            ntot = len(self.freqs[nlist[i]-1])-1
            tauvals = -1.0*(ntot - 2.0*np.arange(len(self.freqs[nlist[i]-1])))/ntot
            probs = self.cdfs[nlist[i]-1]
            plt.scatter(tauvals,probs,color=colors[i])
            # now plot gaussian comparison
            var = 2.0*(2.0*nlist[i]+5.0)/(nlist[i]*(nlist[i]-1)*9.0)
            plt.plot(x,norm.cdf(x,0.0,np.sqrt(var)),color=colors[i])
        plt.legend(nlist,loc='best')
        # plt.savefig('cdfs.png')
        plt.show()

    # use cdfs to return pval
    # default to return two tailed pval
    def pval(self,tau,n,two_tailed=True):
        # enforce tau is between -1 and 1
        if tau < -1.0 or tau > 1.0:
            sys.stderr.write(str(tau)+"\n")
            sys.stderr.write("invalid tau\n")
            #print 'invalid tau'
            return None
        # enforce n is less than our precomputed quantities
        if n > self.N:
            #print 'n is too large'
            sys.stderr.write("n is too large/n")
            return None

        # convert tau to value in terms of number of inversions
        ntot = n*(n-1)/2
        inv_score = int(round((ntot - tau * ntot)/2.0))
        # I'm a little worried about the precision of this,
        # but probably not enough to be really worried for reasonable n
        # since we really only need precision to resolve ntot points

        # if two tailed, we're getting a tail from a symmetric dist
        min_inv_score = min(inv_score,ntot-inv_score)

        if two_tailed:
            pval = self.cdfs[n-1][min_inv_score]*2.0
        else:
            # if one tailed return prob of getting that or fewer inversions
            pval = self.cdfs[n-1][inv_score]

        # if inv_score is 0, might have larger than 0.5 prob
        return min(pval,1.0)



if __name__=="__main__":
    parser = __create_parser__()
    args = parser.parse_args()
    main(args)
