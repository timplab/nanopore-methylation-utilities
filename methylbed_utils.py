import math
import re
import os
import csv
from collections import namedtuple
import re
import pysam
import numpy as np
import pandas as pd

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
    def redo_mcall(self,upper,lower) :
        for i,ratio in enumerate(self.ratios) :
            ratio = float(ratio)
            if ratio > upper :
                call = 1
            elif ratio < lower :
                call = 0
            else :
                call = -1
            self.callarray[i,1] = call
        new_dict = dict()
        for i,key in enumerate(sorted(self.calldict.keys()) ):
            new_dict[key] = methylCall(
                    self.calldict[key].pos,
                    self.callarray[i,1],
                    float(self.ratios[i]),
                    self.seqs[i])
        self.calldict = new_dict

# functions
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

# tabix functions 
def tabix(fpath,window) :
    with pysam.TabixFile(fpath,'r') as tabix :
        entries = [ x for x in tabix.fetch(window)]
    return entries

def tabix_mbed(fpath,regs_df) :
    # fpath is methylation bed file (gzipped and tabix indexed)
    # regs_df must be a pandas data frame with columns : chrom,start,end
    # if only one reg : 
    if not isinstance(regs_df, pd.DataFrame) :
        row = regs_df
        coord = make_coord(row[0],row[1],row[2])
        try :
            raw_list = tabix(fpath,coord)
            out_list = [MethRead(x) for x in raw_list]
        except ValueError : 
            out_list = list()
    else  :
        # if multiple regs, each element of the list is the list of reads in a given reg
        out_list = list()
        for idx,row in regs_df.iterrows() :
            coord = make_coord(row[0],row[1],row[2])
            try :
                raw_list = tabix(fpath,coord)
                out_list.append([MethRead(x) for x in raw_list])
            except ValueError : 
                out_list.append([])
    return out_list

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


