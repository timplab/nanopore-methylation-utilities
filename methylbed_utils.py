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

# sniffles entry
class SnifflesEntry :
    def __init__(self,line) :
        self.line=line.strip()
        self.fields=self.line.split("\t")
        (self.chrom,self.pos,self.id,self.ref,
                self.type,self.qual,self.filter,self.infostring,
                self.format,self.genotype) = self.fields[0:10]
        self.pos = int(self.pos)
        self.type = self.type.strip("<").strip(">")
        self.activate()
    def activate(self) :
        self.parseinfo()
        self.checkcontig()
        self.parsegenotype()
    def checkcontig(self) :
        if "chr" not in self.chrom :
            self.chrom = "chr" + self.chrom
            self.info["CHR2"] = "chr" + self.info["CHR2"]
    def parseinfo(self) :
        self.infofields = [ x.split("=") for x in self.infostring.strip().split(";")]
        self.info = dict()
        self.info["CONFIDENCE"] = self.infofields[0][0]
        for entry in self.infofields[1:] :
            self.info[entry[0]] = entry[1]
        self.info["END"] = int(self.info["END"])
        self.info["RE"] = int(self.info["RE"])
        self.info["SVLEN"] = int(self.info["SVLEN"])
        self.rnames = self.info["RNAMES"].split(',')
    def parsegenotype(self) :
        self.allele= self.genotype.split(":")[0]
        if self.allele== "0/1" :
            self.zygosity = "het"
        elif self.allele == "1/1" :
            self.zygosity = "hom"
        else : self.zygosity = "none"
        self.num_against = int(self.genotype.split(":")[1])
        self.num_for = int(self.genotype.split(":")[2])
        self.coverage = self.num_against+self.num_for

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
    fields=re.split(":|-",coord)
    return fields[0],int(fields[1]),int(fields[2])


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

