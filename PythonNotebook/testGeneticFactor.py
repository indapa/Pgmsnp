#!/usr/bin/env python

from optparse import OptionParser
import pysam
from BedFile import Bedfile
from PileupFactory import PileupFactory
from PgmsnpCommon import *
from common import *
#import pdb
import bx.seq.twobit
from FactorOperations import *
from PgmNetworkFactory import *

""" Let's test contructing our genetic network for the pgmsnp caller
    For simplicity we consider each position in a genome indpendently from each other.
    We iterate through the bedfile interval, and for each position in the interval we will:

    1. construct the network
    2. initialize factors
    3. compute inference to get MAP configuration/probabilities

"""

def main():
    usage = "usage: %prog [options] my.bam "
    parser = OptionParser(usage)
    parser.add_option("--bed", type="string", dest="bedfile", help="bed file with coordinates")
    parser.add_option("--tbf", type="string", dest="tbfile", help=" *.2bit file of the reference sequence")
    parser.add_option("--ped", type="string", dest="pedfile", help= " pedfile of the samples you are analyzing")
    (options, args)=parser.parse_args()
    bamfile=args[0]
    
   
    twobit=bx.seq.twobit.TwoBitFile( open( options.tbfile  ) )
    bedobj=Bedfile(options.bedfile)
    pybamfile=pysam.Samfile( bamfile, "rb" )
    
    # Pfactor gives us a pileup iterator
    Pfactory=PileupFactory(pybamfile,bedobj)




    for pileup_data_obj in Pfactory.yieldPileupData():


        pileup_column_pos=pileup_data_obj.getPileupColumnPos()
        pileup_column_chrom=pileup_data_obj.getChrom()
        refsequence=twobit[pileup_column_chrom][pileup_column_pos:pileup_column_pos+1] #this is the reference base
        print 'Pileup position: ', pileup_data_obj, "\t".join( [pileup_column_chrom,  str(pileup_column_pos), refsequence])
        print

        #lets make our genetic network here:
        pgmNetwork=PgmNetworkFactory(options.pedfile,pileup_column_chrom,pileup_column_pos, refsequence)
        totalSize=pgmNetwork.returnTotalSize()
        """ At some point, we need to construct the clique tree from list of factors.
            But we fight that battle another day ..."""
        #print "pgmNetwork factor list: "
        pgmNetwork.printFactorList()
        pgmFactorList=pgmNetwork.getFactorList()
        print "++++"
        
        for (sample, pileup_sample_data) in pileup_data_obj.yieldSamplePileupData():
            print sample
            sample_idx=pgmNetwork.getSampleNamePedIndex(sample)
            print pgmFactorList[sample_idx + totalSize]
            print
            
            print pileup_sample_data
        
            value=calculateGLL(pileup_sample_data)
            pgmFactorList[sample_idx + totalSize].setVal(value)
            print pgmFactorList[sample_idx + totalSize]
            
            #GLFactor = Factor( [readVar, genoVar], [1,10], [], 'read_phenotype | genotype ')
            #gPrior=LogFactor( returnGenotypePriorFounderFactor(sequence,['A','C','G','T'], genoVar) )
            print
     
        print "==="
        counter=0

if __name__ == "__main__":
    main()
