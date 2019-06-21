#! /usr/bin/env python3 
import os
import sys
import math
import bisect
import argparse
import gzip
import numpy as np
from collections import namedtuple
from methylbed_utils import MethRead,make_coord,bed_to_coord,coord_to_bed
import pysam
from Bio import SeqIO
import re
import multiprocessing as mp
import time
start_time = time.time()

def parseArgs() :
    # dir of source code
    srcpath=sys.argv[0]
    srcdir=os.path.dirname(os.path.abspath(srcpath))
    # parser
    parser = argparse.ArgumentParser(
            description='extract bed entries based on haplotype assignments in a bam file')
    parser.add_argument('-t','--threads',type=int,required=False,default=2, 
            help="number of parallel processes (default : 2 )")
    parser.add_argument('-v','--verbose', action='store_true',default=False,
            help="verbose output")
    parser.add_argument('-b','--bam',type=os.path.abspath,required=True,
            help="bam file - sorted and indexed")
    parser.add_argument('-d','--bed',type=os.path.abspath,required=True,
            help="bed file - sorted, bgzipped, and indexed")
    parser.add_argument('-w','--window',type=str,required=False, 
            help="window from index file to extract [chrom:start-end]")
    parser.add_argument('-r','--regions',type=argparse.FileType('r'),required=False, 
            help="windows in bed file (default: stdin)")
    parser.add_argument('-p','--pre',type=str,required=False,
            default = None,help="output haplotyped bed prefix (default input bam - extension)")
    # parse args
    args = parser.parse_args()
    if not args.pre :
        args.pre = args.bam.rstrip(".bam")
    args.srcdir=srcdir
    return args

# https://stackoverflow.com/questions/107705/disable-output-buffering
class Unbuffered(object):
   def __init__(self, stream):
       self.stream = stream
   def write(self, data):
       self.stream.write(data)
       self.stream.flush()
   def writelines(self, datas):
       self.stream.writelines(datas)
       self.stream.flush()
   def __getattr__(self, attr):
       return getattr(self.stream, attr)
sys.stderr = Unbuffered(sys.stderr)

# https://stackoverflow.com/questions/13446445/python-multiprocessing-safely-writing-to-a-file
def listener(q,prefix,verbose=False) :
    '''listens for messages on the q, writes to file. '''
    if verbose : 
        print("writing output to {}.hap1.bed and {}.hap2.bed".format(prefix,prefix),file=sys.stderr)
    hap1_fp = prefix + ".hap1.bed"
    hap2_fp = prefix + ".hap2.bed"
    hap1_fh = open(hap1_fp,'w')
    hap2_fh = open(hap2_fp,'w')
    qname_list= list()
    def index(a,x) :
        i = bisect.bisect_left(a,x)
        if i != len(a) and a[i] == x :
            return True
        return False
    while True:
        m = q.get()
        if m == 'kill':
            break
        read,hap = m
        fields = read.split("\t")
        qname = ':'.join([fields[0],fields[1],fields[2],fields[3]])
        if not index(qname_list,qname) : 
            bisect.insort(qname_list,qname)
#            if verbose : print("printing {}".format(qname),file=sys.stderr)
            if hap == 1 :
                print(read,file=hap1_fh)
            elif hap == 2 :
                print(read,file=hap2_fh)
#        else :
#            if verbose : print("{} already printed".format(qname),file=sys.stderr)
        q.task_done()
    if verbose : print("total {} reads processed".format(len(qname_list)),file=sys.stderr)
    hap1_fh.close()
    hap2_fh.close()
    q.task_done()

def get_windows_from_bam(bampath,winsize = 100000) :
    with pysam.AlignmentFile(bampath,"rb") as bam :
        stats = bam.get_index_statistics()
        lengths = bam.lengths
        coords = list()
        for stat,contigsize in zip(stats,lengths) : 
            # only fetch from contigs that have mapped reads
            if stat.mapped == 0 :
                continue
            numbins = math.ceil(contigsize/winsize) 
            wins = [ [x*winsize+1,(x+1)*winsize] for x in range(numbins) ] 
            wins[-1][1] = contigsize 
            for win in wins : 
                coords.append(make_coord(stat.contig,win[0],win[1]))
    return coords

def read_tabix(fpath,window) :
    with pysam.TabixFile(fpath) as tabix :
        entries = [ x for x in tabix.fetch(window)]
    rdict = dict()
    # for split-reads, multiple entries are recorded per read name
    for entry in entries:
        qname = entry.split("\t")[3]
        if qname in rdict.keys() :
            rdict[qname].append(entry)
        else : 
            rdict[qname] = [entry]
    return rdict

def split_bed_haplotypes(bampath,bedpath,coord,verbose,q) :
    with pysam.AlignmentFile(bampath,'rb') as bam :
        reads = [x for x in bam.fetch(region=coord)] 
    hap1_qnames = list()
    hap2_qnames = list()
    for read in reads :
        tags = read.tags
        taglabels = [ x[0] for x in tags ]
        if "HP" not in taglabels : continue
        qname = read.query_name
        hap = tags[taglabels.index("HP")][1]
        if hap == 1 :
            hap1_qnames.append(qname)
        elif hap == 2 :
            hap2_qnames.append(qname)
    if len(hap1_qnames) + len(hap2_qnames) == 0 : return
    m_dict = read_tabix(bedpath,coord)
    haplist = [ hap1_qnames,hap2_qnames ]
    num = 0 
    for i in range(len(haplist)) :
        hap = i + 1
        for qname in haplist[i] :
            try : 
                meth_list = m_dict[qname]
            except KeyError : 
                continue
            num += 1
            for meth in meth_list :
                q.put((meth,hap))
    if verbose : 
        print("processed {} bed entries".format(num),file = sys.stderr)

def main() :
    args=parseArgs()
    if args.window : 
        windows = [args.window]
    elif args.regions : 
        windows = [ bed_to_coord(x) for x in args.regions ]
    else :
        if args.verbose : 
            print("extracting all reads in the bam file",file=sys.stderr)
        windows = get_windows_from_bam(args.bam,100000)
    if args.verbose : print("{} regions to parse".format(len(windows)),file=sys.stderr)
    # initialize mp
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(processes=args.threads)
    if args.verbose : print("using {} parallel processes".format(args.threads),file=sys.stderr)
    # watcher for output
    watcher = pool.apply_async(listener,(q,args.pre,args.verbose))
    # start processing
    jobs = [ pool.apply_async(split_bed_haplotypes,
        args = (args.bam,args.bed,win,args.verbose,q))
        for win in windows ]
    output = [ p.get() for p in jobs ]
    # done
    q.put('kill')
    q.join()
    pool.close()
    if args.verbose : print("time elapsed : {} seconds".format(time.time()-start_time),file=sys.stderr)

if __name__=="__main__":
    main()

