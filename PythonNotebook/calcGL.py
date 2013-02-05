#!/usr/bin/env python
import os
import string
import re
from optparse import OptionParser
import pysam
import sys
from BedFile import Bedfile
from PileupData import PileupData
from PileupFactory import PileupFactory
from PgmsnpCommon import *
from common import *
import pdb
import bx.seq.twobit
from FactorOperations import *

def main():
    usage = "usage: %prog [options] my.bam "
    parser = OptionParser(usage)
    parser.add_option("--bed", type="string", dest="bedfile", help="bed file with coordinates")
    parser.add_option("--tbf", type="string", dest="tbfile", help=" *.2bit file of the reference sequence")
    (options, args)=parser.parse_args()
    bamfile=args[0]
    
    
    twobit=bx.seq.twobit.TwoBitFile( open( options.tbfile  ) )
    bedobj=Bedfile(options.bedfile)
    pybamfile=pysam.Samfile( bamfile, "rb" )
    #pyfai=pysam.Fastafile(options.faidxfile)
    Pfactory=PileupFactory(pybamfile,bedobj)
    counter=0
    for pileup_data_obj in Pfactory.yieldPileupData():


        pileup_column_pos=pileup_data_obj.getPileupColumnPos()
        pileup_column_chrom=pileup_data_obj.getChrom()
        sequence=twobit[pileup_column_chrom][pileup_column_pos:pileup_column_pos+1]
        print 'Pileup position: ', pileup_data_obj, "\t".join( [pileup_column_chrom,  str(pileup_column_pos), sequence])
        print
        for (sample, pileup_sample_data) in pileup_data_obj.yieldSamplePileupData():
            #print pileup_sample_data
            genoVar=counter
            counter+=1
            readVar=counter
            counter+=1
            value=calculateGLL(pileup_sample_data)
            GLFactor = Factor( [readVar, genoVar], [1,10], [], 'read_phenotype | genotype ')
            gPrior=LogFactor( returnGenotypePriorFounderFactor(sequence,['A','C','G','T'], genoVar) )
            
            #print sample, value
            GLFactor.setVal(value)
            print
            print gPrior
            print
            print GLFactor
            print 
            print FactorSum(gPrior, GLFactor)
            print 
        print "==="
        counter=0

if __name__ == "__main__":
    main()
