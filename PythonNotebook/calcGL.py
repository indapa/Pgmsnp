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

def main():
    usage = "usage: %prog [options] my.bam "
    parser = OptionParser(usage)
    parser.add_option("--bed", type="string", dest="bedfile", help="bed file with coordinates")

    (options, args)=parser.parse_args()
    bamfile=args[0]
    bedobj=Bedfile(options.bedfile)
    pybamfile=pysam.Samfile( bamfile, "rb" )
    Pfactory=PileupFactory(pybamfile,bedobj)
    counter=0
    for pileup_obj in Pfactory.yieldPileupData():

        print 'Pileup position: ', pileup_obj
        for (sample, pileup_sample_data) in pileup_obj.yieldSamplePileupData():
            genoVar=counter
            counter+=1
            readVar=counter
            counter+=1
            value=calculateGLL(pileup_sample_data)
            GLFactor = Factor( [readVar, genoVar], [1,10], [], 'read_phenotype | genotype ')
            #print sample, value
            GLFactor.setVal(value)
            print GLFactor
            print 
        print "==="
        counter=0

if __name__ == "__main__":
    main()
