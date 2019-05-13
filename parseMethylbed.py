#! /usr/bin/env python
import os
import math
import bisect
import sys
import argparse
import gzip
import numpy as np
from collections import namedtuple
from methylbed_utils import MethRead

def parseArgs() :
    # dir of source code
    srcpath=sys.argv[0]
    srcdir=os.path.dirname(os.path.abspath(srcpath))
    # parser
    parser = argparse.ArgumentParser(description='parse methylation bed files')
    subparsers = parser.add_subparsers()
    # parent parser
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-i','--input',type=os.path.abspath,required=False,
            help="input methylation bed file (default stdin)")
    parent_parser.add_argument('-o','--out',type=argparse.FileType('w'),required=False,
            default=sys.stdout,help="output file path (default stdout)")
    parent_parser.add_argument('-v','--verbose',action='store_true',default=False,
            help="verbose output")
    parent_parser.add_argument('-m','--mod',type=str,required=False,
            help="methylation motif (one of 'cpg' and 'gpc', default cpg)",default="cpg")
    # parser for frequency
    parser_freq = subparsers.add_parser('frequency',parents=[parent_parser],
            help = 'get frequency')
    parser_freq.set_defaults(func=getFreq)
    # parser for readlevel methylation
    parser_readlevel = subparsers.add_parser('readlevel',parents=[parent_parser],
            help = 'get read-level methylation')
    parser_readlevel.set_defaults(func=getReadlevel)
    # parser for intersect methylation
    parser_intersect = subparsers.add_parser('intersect',parents=[parent_parser],
            help = 'get intersect methylation')
    parser_intersect.set_defaults(func=getMethIntersect)
    # parser for per-read average methylation
    parser_per_read = subparsers.add_parser('per_read',parents=[parent_parser],
            help = 'get per-read average methylation')
    parser_per_read.set_defaults(func=getPerReadAverage)
    parser_per_read.add_argument('-e','--exclude',type=str,required=False,default="GCG",
            help = 'motif to exclude (default "GCG")')

    # parse args
    args = parser.parse_args()
    args.srcdir=srcdir
    if args.mod=="cpg":
        args.motif="CG"
    elif args.mod=="gpc":
        args.motif="GC"
    return args

def printsite(stats,out):
    print("{}\t{}\t{}\t{}".format(stats.rname,stats.pos,
        stats.num_reads,stats.num_methylated),file=out)
def printcyto(stats,mod,out):
    if mod=="cpg":
        motif="CG"
    elif mod =="gpc" :
        motif="GC"
    conind=stats.seq.index(motif)
    context=stats.seq[conind:conind+4]
    print("{}\t{}\t*\t{}\t{}\tCG\t{}".format(stats.rname,stats.pos,
        stats.num_methylated,
        stats.num_reads-stats.num_methylated,
        context),file=out)

class SiteStats:
    def __init__(self,methcall,chrom):
        self.rname=chrom
        self.pos=methcall.pos
        self.num_reads = 0
        self.num_methylated = 0
        self.seq=methcall.seq
    def update(self,methcall):
        if methcall.call==-1:
            return
        self.num_reads+=1
        self.num_methylated+=methcall.call
    def printFreq(self,context,out):
        if self.num_reads == 0 : return # all calls here were bad
        if context == "CG" :
            motifidx=self.seq.index(context,1)
            tricontext=self.seq[motifidx-1:motifidx+2]
        elif context == "GC" :
            motifidx=self.seq.index(context)
            tricontext=self.seq[motifidx:motifidx+3]
        motifcontext=self.seq[motifidx:motifidx+2]
        print("\t".join([str(x) for x in [
            self.rname,
            self.pos+1,
            "*",
            self.num_methylated,
            self.num_reads-self.num_methylated,
            motifcontext,tricontext]]),file=out)

def getFreq(args,in_fh):
    if args.verbose : print("getting frequency",file=sys.stderr)
    sites=dict()
    n = 0
    for line in in_fh:
        try : 
            line = line.decode('ascii')
        except AttributeError :
            pass
        n += 1
        read=MethRead(line)
#        # debug
#        print(line,file=sys.stdout)
#        print(read.ratios)
#        print(read.keys)
#        #
        sitekeys=sorted(sites.keys())
#        print(sitekeys)
        try : 
            # print everything in sites if chromosome is different
            if read.rname != sites[sitekeys[0]].rname :
                printind=len(sitekeys)
            else : 
                # get index of new position
                printind=bisect.bisect_left(sitekeys,read.keys[0])
        except IndexError : 
            printind=0
        if printind != 0 :
            for i in range(printind):
                key=sitekeys[i]
                sites[key].printFreq(args.motif,args.out)
                sites.pop(key)
        for key in read.keys:
            if  key not in sites.keys():
                sites[key] = SiteStats(read.calldict[key],read.rname)
            sites[key].update(read.calldict[key])
        if args.verbose : 
            if n % 10000  == 0:
                print("parsed {} lines".format(n),file=sys.stderr)
    for key in sorted(sites.keys()):
        sites[key].printFreq(args.motif,args.out)

def getReadlevel(args,in_fh):
    for line in in_fh:
        try : 
            line = line.decode('ascii')
        except AttributeError :
            pass
        read=MethRead(line)
        callkeys=sorted(read.calldict)
        for key in callkeys:
            methcall=read.calldict[key]
            outlist=[read.rname]+[str(x) for x in methcall]
            outlist.insert(3,read.qname)
            print("\t".join(outlist),file=args.out)

def getMethIntersect(args,in_fh) :
    for line in in_fh : 
        try : 
            line = line.decode('ascii')
        except AttributeError :
            pass
        intersect=methInt(line)
        read=intersect.methread
        region=intersect.region
        start=region.start
        end=region.end
        for key in read.keys :
            methcall=read.calldict[key]
            if key >= start and key < end :
                outlist=([region.rname]+
                        [str(x) for x in methcall]+
                        [str(x) for x in region.items]+
                        [str(methcall.pos-region.start)])
                outlist.insert(2,str(methcall.pos+1))
                outlist.insert(4,read.qname)
                print("\t".join(outlist),file=args.out)

def getPerReadAverage(args,in_fh) :
    for line in in_fh :
        try : 
            line = line.decode('ascii')
        except AttributeError :
            pass
        read = MethRead(line)
        seqs = [x.seq for x in read.calldict.values()]
        # filter out calls with unwanted seq
        filteridx = [idx for idx,s in enumerate(seqs) if args.exclude not in s]
        calls_all = np.array([x.call for x in read.calldict.values()])[filteridx]
        ratios = np.array([x.ratio for x in read.calldict.values()])[filteridx]
        # filter out unconfident calls
        calls = [ x for x in calls_all if x != -1 ]
        # drop empty call strings
        if len(calls) > 0 :
            print('\t'.join([str(x) for x in [
                read.rname,
                read.start,
                read.end,
                read.qname,
                np.mean(ratios),'.',
                np.mean(calls),
                len(calls)]]),file=args.out)
      
def main() :
    args=parseArgs()
    # read in input
    if args.input:
        if args.input.split('.')[-1]=="gz":
            in_fh=gzip.open(args.input,'rt')
        else :
            in_fh = open(args.input,'r')
    else:
        in_fh = sys.stdin
    args.func(args,in_fh)
    in_fh.close()

if __name__=="__main__":
    main()

