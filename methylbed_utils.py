import math
import re
import os
import csv
from collections import namedtuple
import re
import pysam
import numpy as np

# for read from a methylation bed file
methylCall = namedtuple('methylCall', ['pos','call','ratio','seq'])
class MethRead :
    def __init__(self,string):
        self.string=string
        self.fields=string.strip().split("\t")
        self.rname=self.fields[0]
        self.start=int(self.fields[1])
        self.end=int(self.fields[2])
        self.rlen=self.end-self.start
        self.qname=self.fields[3]
        self.methstring=self.fields[4]
        if self.fields[5] == "+" or self.fields[5] == "-" :
            self.strand = self.fields[5]
            self.fields[5] = self.fields[6]
            self.fields[6] = self.fields[7]
        self.ratios=self.fields[5].strip().split(",")
        self.seqs=self.fields[6].strip().split(",")
        self.calldict=self.parseMeth()
        self.keys=sorted(self.calldict.keys())
        self.callarray=self.getArray(self.calldict)
    def make_key(self,pos):
        return int(pos)
    def parseMeth(self):
        calldict=dict()
        pos=int(self.start)
        calls=re.findall('(\d+)([umx])?',self.methstring)
#        # temporary fix for when first distance is not 0
#        if int(calls[0][0]) != 0 :
#            calls = [('0','x')] + calls
#            self.ratios = ['0'] + self.ratios
#            self.seqs = ['GCG'] + self.seqs
#        ##
        for i,(dist,call) in enumerate(calls):
            pos+=int(dist)
            if call=="x" :
                is_meth=-1
            else :
                is_meth=int(call=="m")
            calldict[self.make_key(pos)]=methylCall(
                    pos,
                    is_meth,
                    float(self.ratios[i]),
                    self.seqs[i])
        return calldict
    def getArray(self,calldict) :
        callarray=np.array([(x,calldict[x].call) for x in sorted(calldict.keys())])
        return callarray

# functions 
def tabix(fpath,window) :
    with pysam.TabixFile(fpath,'r') as tabix :
        entries = [ x for x in tabix.fetch(window)]
    return entries
    
def make_coord(chrom,start,end) :
    if start < 1 : start = 1
    return chrom+":"+str(start)+"-"+str(end)

def bed_to_coord(bedentry) :
    fields=bedentry.strip().split("\t")
    start = str(int(fields[1])+1)
    return fields[0]+":"+start+"-"+fields[2]

def coord_to_bed(coord) :
    fields = coord.split(":")
    numbers = fields[1].split("-")
    return fields[0],int(numbers[0]),int(numbers[1])


def read_bam(fpath,window) :
    with pysam.AlignmentFile(fpath,'rb') as bam :
        bam_entries = [ x for x in bam.fetch(region=window) ]
    bamdict = dict()
    for bam in bam_entries :
        try :
            bamdict[bam.query_name].append(bam)
        except :
            bamdict[bam.query_name] = [bam]
    return bamdict

