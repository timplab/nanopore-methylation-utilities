#! /usr/bin/env python3
import sys
import gzip
import csv
import argparse
import re
import pysam

def parseArgs():
    parser = argparse.ArgumentParser( description='Generate methylation bedGraph file')
    parser.add_argument('-v', '--verbose', action="store_true",required=False, default = False,
            help="verbose output")
    parser.add_argument('-c', '--call-threshold', type=float, required=False, default=1.5,
            help="absolute value of threshold for methylation call (default : 1.5)")
    parser.add_argument('-i', '--input', type=str, required=False,help="input methylation tsv file (default stdin)")
    parser.add_argument('-q', '--mod',type=str,required=False,default='cpg',help="modification motif; one of cpg,gpc,dam,cpggpc")
    parser.add_argument('-g','--genome',type=str,required=True,help="Reference genome fasta")
    parser.add_argument('-e', '--exclude',type=str,required=False,help="motif to exclude from reporting")
    parser.add_argument('-w', '--window',type=int,required=False,default=2,
            help="number of nucleotides to report on either side of called nucleotide")
    parser.add_argument('--nome',action="store_true",required=False,default=False,
            help="nanoNOMe - remove calls at GCG motifs")
    parser.add_argument('--offset',type=int,required=False,default=0,
            help="nanopolish coordinate offset (1-based)")
    args = parser.parse_args()
    assert(args.call_threshold is not None)
    return args

class readQuery:
    def __init__(self,record,thr,motif,offset,win_size,nome):
        self.thr = thr
        self.motif = motif
        self.offset = offset
        self.win_size = win_size
        self.nome = nome
        self.qname = record['read_name']
        self.rname = record['chromosome']
        self.start = int(record['start']) + self.offset
        self.end = self.start + 1
        self.dist=[]
        self.seq=[]
        self.ratio=[]
        self.call=[]

    def update(self,record,seq_dict):
        llr=float(record['log_lik_ratio'])
        call=self.call_methylation(llr) # call methylation
#        start = int(record['start']) + self.offset
#        sequence = record['sequence']
        start = int(record['start'])-self.win_size
        sequence = seq_dict[self.rname][start:int(record['end'])+2+self.win_size].upper()
        pos = sequence.find(self.motif)
        pos_list = []
        seq_list = []
        # read groups
        while pos != -1 :
            seq = sequence[pos-self.win_size:pos+self.win_size+1]
            if not self.nome or ( "GCG" not in seq ):
                pos_list.append(pos)
                self.ratio.append(llr)
                self.call.append(call)
                seq_list.append(seq)
            pos = sequence.find(self.motif,pos+1)
        if len(pos_list) == 0 : return
        coord_list = [ x+start for x in pos_list ]
        # get distance (CIGAR style)
        dist_list = [ x-y for x,y in 
                zip(coord_list,[self.end-1]+coord_list) ]
        # update 
        self.end = coord_list[-1] + 1
        self.seq = self.seq + seq_list
        self.dist = self.dist + dist_list
    def call_methylation(self,llr):
        if abs(llr) < self.thr :
            call="x"
        elif llr > 0 :
            call="m"
        else : 
            call="u"
        return call
    def summarizeCall(self):
        summary=""
        for ind,call in zip(self.dist,self.call):
            summary=summary+"{}{}".format(ind,call)
        return summary
    def printRead(self,extra=""):
        if len(self.call) == 0 : return # after filtering GCG there is no data in this read
        # if GCG is first motif, adjust start to make first index 0
        if self.dist[0] != 0 :
            self.start += self.dist[0]
            self.dist[0] = 0
        # when motif is GC, positions need to shift + one
        if self.motif == "GC" :
            self.start += 1
            self.end += 1 
        def catList(strlist,delim=""):
            return delim.join([str(x) for x in strlist])
        print("\t".join([self.rname,str(self.start),str(self.end),
            self.qname,
            self.summarizeCall(),
            catList(self.ratio,","),
            catList(self.seq,",")])+extra)

def summarizeMeth(args):
    # first load ref into memory to speed up lookup
    if args.verbose : print("loading reference sequence",file = sys.stderr)
    fa = pysam.FastaFile(args.genome)
    contigs = fa.references
    seq_dict = { x:fa.fetch(x) for x in contigs }
    # motif and offset based on modification
    if args.mod=="cpg" :
        motif="CG"
    elif args.mod=="gpc" :
        args.offset += 1  # coordinate of C is +1 since G comes before C
        motif="GC"
    elif args.mod=="dam" :
        args.offset += 1  # coordinate of C is +1 since G comes before C
        motif="GATC"
        args.window += 2
    # read in data
    if args.input:
        if args.input.strip().split(".")[-1] == "gz" :
            in_fh = gzip.open(args.input,"rt")
        else : 
            in_fh = open(args.input)
    else:
        in_fh = sys.stdin
    csv_reader = csv.DictReader(in_fh,delimiter='\t')
    if args.verbose : print("processing calls",file = sys.stderr)
    if args.mod != "cpggpc" :
        for record in csv_reader:
            # skip queries that have an undesired motif in the call group
            if args.exclude:
                if args.exclude in record['sequence'] :
                    continue
            # update query
            try :
                if ((record['read_name'] != read.qname) or
                        (record['chromosome'] != read.rname) or
                        (int(record['start']) < read.end)):
                    read.printRead()
                    read = readQuery(record,args.call_threshold,
                            motif,args.offset,args.window,args.nome) 
            except NameError : # has not been initialized
                read = readQuery(record,args.call_threshold,
                        motif,args.offset,args.window,args.nome)
            except ValueError : # header or otherwise somehow faulty
                continue
            read.update(record,seq_dict)
        # finishing up - print the unprinted read
        if read.qname : 
            read.printRead()
    else :
        # if motif is cpggpc, the output is split between cpg calls and gpc calls
        read = dict()
        for record in csv_reader :
            if record['chromosome'] not in contigs :
                if args.verbose : 
                    print("{} not in reference, skipping.".format(record['chromosome']),file = sys.stderr)
                continue
            motif = record['motif']
            if record['read_name'] == "read_name" : continue # skip extra headers
            record['sequence'] = re.sub("M","C",record['sequence'])
            # update query
            try :
                if ((record['read_name'] != read[motif].qname) or
                        (record['chromosome'] != read[motif].rname) or
                        (int(record['start']) < read[motif].end)):
                    for key in read.keys() :
                        read[key].printRead("\t"+key)
                    read = dict()
                    read[motif] = readQuery(record,args.call_threshold,
                            motif,args.offset,args.window,args.nome)
            except KeyError : # has not been initialized
                read[motif] = readQuery(record,args.call_threshold,
                        motif,args.offset,args.window,args.nome)
            except ValueError : # header or otherwise somehow faulty
                continue
            read[motif].update(record,seq_dict)
        # finishing up - print the unprinted read
        for key in read.keys() :
            read[key].printRead("\t"+key)

    in_fh.close()

if __name__=="__main__":
    args=parseArgs()
    summarizeMeth(args)
