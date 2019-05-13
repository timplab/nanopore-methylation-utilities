#! /usr/bin/env python
import sys
import csv
import argparse

def parseArgs():
    parser = argparse.ArgumentParser( description='Generate methylation bedGraph file')
    parser.add_argument('-c', '--call-threshold', type=float, required=False, default=2.5,
            help="absolute value of threshold for methylation call (default : 2.5)")
    parser.add_argument('-i', '--input', type=str, required=False,help="input methylation tsv file (default stdin)")
    parser.add_argument('-m', '--mod',type=str,required=False,default='cpg',help="modification motif; one of cpg,gpc,dam")
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
        self.strand = record['strand']
        self.start = int(record['start']) + self.offset
        self.end = self.start + 1
        self.dist=[]
        self.seq=[]
        self.ratio=[]
        self.call=[]
        self.update(record)

    def update(self,record):
        llr=float(record['log_lik_ratio'])
        call=self.call_methylation(llr) # call methylation
        start = int(record['start']) + self.offset
        sequence = record['sequence']
        pos = sequence.find(self.motif)
        pos_list = []
        seq_list = []
        # read groups
        while pos != -1 :
            seq = sequence[pos-self.win_size:pos+self.win_size+1]
            if not self.nome or "GCG" not in seq:
                pos_list.append(pos)
                self.ratio.append(llr)
                self.call.append(call)
                seq_list.append(seq)
            pos = sequence.find(self.motif,pos+1)
        if len(pos_list) == 0 : return
        coord_list = [ x+start-pos_list[0] for x in pos_list ]
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
    def printRead(self):
        if len(self.call) == 0 : return # after filtering GCG there is no data in this read
        # if GCG is first motif, adjust start to make first index 0
        if self.dist[0] != 0 :
            self.start += self.dist[0]
            self.dist[0] = 0
        def catList(strlist,delim=""):
            return delim.join([str(x) for x in strlist])
        print("\t".join([self.rname,str(self.start),str(self.end),
            self.qname,
            self.summarizeCall(),
            catList(self.ratio,","),
            catList(self.seq,",")]))

def summarizeMeth(args):
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
        in_fh = open(args.input)
    else:
        in_fh = sys.stdin
    csv_reader = csv.DictReader(in_fh,delimiter='\t')
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
            else :
                read.update(record)
        except NameError : # has not been initialized
            read = readQuery(record,args.call_threshold,
                    motif,args.offset,args.window,args.nome)
        except ValueError : # header or otherwise somehow faulty
            continue
    # finishing up - print the unprinted read
    if read.qname : 
        read.printRead()
    in_fh.close()

if __name__=="__main__":
    args=parseArgs()
    summarizeMeth(args)
