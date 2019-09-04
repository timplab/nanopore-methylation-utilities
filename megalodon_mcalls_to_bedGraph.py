#! /usr/bin/env python
import sys
import gzip
import csv
import argparse
import pysam

def parseArgs():
    parser = argparse.ArgumentParser( description='Generate methylation bedGraph file')
    parser.add_argument('-v', '--verbose', action="store_true",required=False, default = False,
            help="verbose output")
    parser.add_argument('-i', '--input', type=str, required=False,help="input methylation tsv file (default stdin)")
    parser.add_argument('-g', '--genome', type=str, required=True,help="reference fasta file, indexed")
    parser.add_argument('-q', '--mod',type=str,required=False,default='cpg',help="modification motif; one of cpg,gpc")
    parser.add_argument('-c', '--call-threshold', type=float, required=False, default=2,
            help="absolute value of threshold for methylation call (default : 2)")
    parser.add_argument('-e', '--exclude',type=str,required=False,help="motif to exclude from reporting")
    parser.add_argument('--nome',action="store_true",required=False,default=True,
            help="nanoNOMe - retain calls at CpG and GpC sites and remove calls at GCG motifs - currently default true")
    parser.add_argument('-w', '--window',type=int,required=False,default=5,
            help="number of nucleotides to report on either side of called nucleotide")
    args = parser.parse_args()
    assert(args.call_threshold is not None)
    return args

# https://codereview.stackexchange.com/questions/151329/reverse-complement-of-a-dna-string
def revcomp(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in dna[::-1]])

class readQuery:
    def __init__(self,record,thr,motif) :
        self.thr = thr
        self.motif = motif
        self.qname = record['read_id']
        self.rname = record['chrm']
        self.strand = record['strand']
        if self.strand == 1 :
            self.update = self.update_positive
        else :
            self.update = self.update_negative
        self.start = record['pos']
        self.end = self.start + 1
        self.dist=[]
        self.seq=[]
        self.ratio=[]
        self.call=[]
        self.update(record)

    def update_positive(self,record):
        ratio=round(float(record['mod_log_prob']) - float(record['can_log_prob']),5)
        call=self.call_methylation(ratio) # call methylation
        pos = record['pos']
        sequence = record['sequence']
        # get distance (CIGAR style)
        dist = pos - ( self.end - 1 )
        # update 
        self.end = pos + 1
        self.call.append(call)
        self.ratio.append(ratio)
        self.seq.append(sequence)
        self.dist.append(dist)
    def update_negative(self,record):
        ratio=round(float(record['mod_log_prob']) - float(record['can_log_prob']),5)
        call=self.call_methylation(ratio) # call methylation
        pos = record['pos']
        sequence = record['sequence']
        # get distance (CIGAR style)
        dist = self.start - pos
        # update 
        self.start = pos
        self.call = [ call ] + self.call
        self.ratio= [ ratio ] + self.ratio
        self.seq = [ sequence ] + self.seq
        self.dist = [ dist ] + self.dist

    def call_methylation(self,ratio):
        if abs(ratio) < self.thr :
            call="x"
        elif ratio > 0 :
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
        # get strand
        if self.strand == 1 :
            strand = "+"
        else : 
            strand = "-"
            self.dist = [ 0 ] + self.dist[:-1] # negative strand special case - cigar order is offset
        # if GCG is first motif, adjust start to make first index 0
        if self.dist[0] != 0 :
            self.start += self.dist[0]
            self.dist[0] = 0
        def catList(strlist,delim=""):
            return delim.join([str(x) for x in strlist])
        print("\t".join([self.rname,str(self.start),str(self.end),
            self.qname,
            self.summarizeCall(),
            strand,
            catList(self.ratio,","),
            catList(self.seq,",")]))

def summarizeMeth(args):
    # motif 
    if args.mod == "cpg" :
        motif = "CG"
    elif args.mod == "gpc" :
        motif = "GC"
    # nanonome options
    if args.nome :
        args.exclude = "GCG"
    # first load ref into memory to speed up lookup
    if args.verbose : print("loading reference sequence",file = sys.stderr)
    fa = pysam.FastaFile(args.genome)
    contigs = fa.references
    seq_dict = { x:fa.fetch(x) for x in contigs }
    # read in data
    if args.input:
        if args.input.strip().split(".")[-1] == "gz" :
            in_fh = gzip.open(args.input,"rt")
        else : 
            in_fh = open(args.input)
    else:
        in_fh = sys.stdin
    if args.verbose : print("parsing data",file = sys.stderr)
    csv_reader = csv.DictReader(in_fh,delimiter='\t')
    for record in csv_reader:
        pos = int(record['pos'])
        record['strand'] = int(record['strand'])
        seq = seq_dict[record['chrm']][pos-1:pos+2].upper()
       	# reverse complement if neg strand
        if record['strand'] == -1 :
            seq = revcomp(seq)
        # disregard non-motif and undesired motifs
        if ((motif not in seq ) or (seq == args.exclude)) :
            continue
        # get seq context
        context = seq_dict[record['chrm']][pos-args.window:pos+args.window+1].upper()
        record['sequence'] = context
        record['pos'] = pos
        # update query
        try :
            if ((record['read_id'] != read.qname) or
                    (record['chrm'] != read.rname) or
                    (record['strand'] != read.strand) or
                    (record['strand'] == 1 and pos < read.end) or
                    (record['strand'] == -1 and pos > read.start)):
                read.printRead()
                read = readQuery(record,args.call_threshold,motif)
            else :
                read.update(record)
        except NameError : # has not been initialized
            read = readQuery(record,args.call_threshold,motif)
        except ValueError : # header or otherwise somehow faulty
            continue
    # finishing up - print the unprinted read
    if read.qname : 
        read.printRead()
    in_fh.close()

if __name__=="__main__":
    args=parseArgs()
    summarizeMeth(args)
