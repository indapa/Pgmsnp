#!/usr/bin/env python
from optparse import OptionParser
import pysam
from BedFile import Bedfile
from PileupFactory import PileupFactory
from PgmsnpCommon import *
from common import *
#import networkx as nx
import bx.seq.twobit
from FactorOperations import *
from PgmNetworkFactory import *
from CliqueTreeOperations import *
from PedigreeFactors import Pedfile
import datetime
from collections import defaultdict
from itertools import chain
#import matplotlib.pyplot as plt
import pdb
import os
import subprocess
""" Let's test contructing our genetic network for the pgmsnp caller
    For simplicity we consider each position in a genome indpendently from each other.
    We iterate through the bedfile interval, and for each position in the interval we will:

    1. construct the network
    2. initialize factors
    3. execute  inference to get MAP configuration/probabilities

"""

def main():
    """ getting the git hash in python: http://stackoverflow.com/a/14989911/1735942"""
    label = subprocess.check_output(["git", "rev-parse", "HEAD"])
    today=datetime.datetime.today()
    datestr=today.strftime("20%y%m%d")
    usage = "usage: %prog [options] my.bam "
    parser = OptionParser(usage)
    parser.add_option("--bed", type="string", dest="bedfile", help="bed file with coordinates")
    parser.add_option("--tbf", type="string", dest="tbfile", help=" *.2bit file of the reference sequence")
    parser.add_option("--ped", type="string", dest="pedfile", help= " pedfile of the samples you are analyzing")
    parser.add_option("--min-nonref-count", dest="minAlt", default=2, help="minimum observation of nonref allele for genotype inference (default 2)")
    (options, args)=parser.parse_args()
    
    commandline=";".join(sys.argv)
    
    bamfile=args[0]
    bamfilebasename=return_file_basename(bamfile)
    vcfoutput=".".join(['pgmsnp',bamfilebasename, datestr, 'vcf'])
    vcfh=open(vcfoutput,'w')
    ALPHABET=['AA', 'AC', 'AG', 'AT', 'CC', 'CG', 'CT', 'GG', 'GT', 'TT']
   
    twobit=bx.seq.twobit.TwoBitFile( open( options.tbfile  ) )
    bedobj=Bedfile(options.bedfile)
    pybamfile=pysam.Samfile( bamfile, "rb" )
    pedfile=Pedfile(options.pedfile)
    pedfile.parsePedfile()

    samplenames=pedfile.returnIndivids()
    samplestring="\t".join(samplenames)
    

    #vcf metainfolines
    vcf_metalines=[]
    vcf_metalines.append ( "##fileformat=VCFv4.1")
    vcf_metalines.append( "##fileDate="+datestr )
    vcf_metalines.append("##testGeneticFactor="+commandline)
    vcf_metalines.append("##version="+label.strip())
    vcf_metalines.append("##reference="+options.tbfile)
    vcf_metalines.append("##pedfile="+options.pedfile)
    vcf_metalines.append("##bedfile="+options.bedfile)
    vcf_metalines.append("##bamfile="+bamfile)
    vcf_metalines.append("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">")
    vcf_metalines.append( "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" )
    vcf_metalines.append( "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" )
    vcf_metalines.append( "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">" )
    vcf_metalines.append( "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" )
    vcf_column_headers=["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",samplestring]
    vcf_metalines.append( "\t".join(vcf_column_headers))

    vcf_metaline_output="\n".join(vcf_metalines)
    vcfh.write(vcf_metaline_output+"\n")
    #print vcf_metaline_output



    #print samplenames

    #return the complete enumeration of all possible offspring genotype priors
    punnetValues=returnPunnetValues(4)

    # Pfactor gives us a pileup iterator
    Pfactory=PileupFactory(pybamfile,bedobj)



    ##CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",samplestring

    for pileup_data_obj in Pfactory.yieldPileupData():
        
        refDepth=defaultdict(lambda: defaultdict(int))
        altDepth=defaultdict(lambda: defaultdict(int))
        pileup_column_pos=pileup_data_obj.getPileupColumnPos()
        pileup_column_chrom=pileup_data_obj.getChrom()
        refsequence=twobit[pileup_column_chrom][pileup_column_pos:pileup_column_pos+1] #this is the reference base
        refDepth=defaultdict(lambda: defaultdict(int))
        altDepth=defaultdict(lambda: defaultdict(int))
        skipsite=False
        #print 'Pileup position: ', pileup_data_obj, "\t".join( [pileup_column_chrom,  str(pileup_column_pos), refsequence])
        #print
        
        qual="."
        filter='.'
        siteDP="DP="+str(pileup_data_obj.getPileupDepth())
        altDP=0
        sampleDepth=defaultdict(int)


        #lets make our genetic network here:
        pgmNetwork=PgmNetworkFactory(options.pedfile,pileup_column_chrom,pileup_column_pos, refsequence,punnetValues)
        totalSize=pgmNetwork.returnTotalSize()
        #print "totalSize: ", totalSize
        observedSamples=[]
        sampleNames=set(pgmNetwork.getSampleNames())
        #print "pgmNetwork factor list: "
        #pgmNetwork.printFactorList()
        pgmFactorList=pgmNetwork.getFactorList()
        #print "++++"
        
        for (sample, pileup_sample_data) in pileup_data_obj.yieldSamplePileupData():
            sampleDepth[sample]=len(pileup_sample_data)
            observedSamples.append(sample)
            sample_idx=pgmNetwork.getSampleNamePedIndex(sample)
            
            for data in pileup_sample_data:
                if data.basecall != refsequence: altDepth[data.sample][data.basecall]+=1
                if data.basecall == refsequence: refDepth[data.sample][data.basecall]+=1
            
       
            
            value=calculateGLL(pileup_sample_data)
            

            pgmFactorList[sample_idx + totalSize].setVal(value)
            #we initially get the log-likelihood, but go back to probablity space first
            pgmFactorList[sample_idx + totalSize] = ExpFactor( pgmFactorList[sample_idx + totalSize] )
            #print pgmFactorList[sample_idx + totalSize]
            
            #GLFactor = Factor( [readVar, genoVar], [1,10], [], 'read_phenotype | genotype ')
            #gPrior=LogFactor( returnGenotypePriorFounderFactor(sequence,['A','C','G','T'], genoVar) )
            #print
        
    
        if sum(chain.from_iterable(d.itervalues() for d in altDepth.itervalues())) < options.minAlt: 
       #print pileup_data_obj.getPileupColumnPos(), altDP
            continue
        
        observedSamples=set(observedSamples)
        unobservedSamples=sampleNames-observedSamples
        unobservedIdxs=[ pgmNetwork.getSampleNamePedIndex(sample) for sample in unobservedSamples  ]
        #for samples lacking read data, delete them from the list of
        for idx in unobservedIdxs:
            #del pgmFactorList[idx + totalSize]
            value=calculateNoObservationsGL()
            #pdb.set_trace()
            pgmFactorList[idx + totalSize].setVal(value)

        #print pgmFactorList
        #for f in pgmFactorList:
        #    print f
        #continue
        #print "====="
        siteNS="NS="+str(len(observedSamples))
        infoString=";".join([siteNS, siteDP])
        
        
        
        #for f in pgmFactorList:
        #    print f
        #    print

        #pdb.set_trace()
        #cTree=CreatePrunedInitCtree(pgmFactorList)
        #print 'cTree constructed, pruned, and initalized potential'


        #G=nx.from_numpy_matrix( cTree.getEdges() )
        #nx.draw_shell(G)
        #plt.show()


        #print cTree

        #get the max marginal factors
        MAXMARGINALS=ComputeExactMarginalsBP(pgmFactorList, [], 1)
        
        #MARGINALS= ComputeExactMarginalsBP( pgmFactorList)
        #this log normalizes the data
        for m in MAXMARGINALS:
            #print m.getVal()
            m.setVal( np.log( lognormalize(m.getVal()   )   ) )
            #print np.sum( lognormalize(m.getVal() ) )
            #print'lognormalize: ',  lognormalize(m.getVal() )
            #print m.getVal()
            #print
        #    print  m.getVal()
        #    print
        #print "==="
        #get the max value of each factor
        
        
                
                           
        MAXvalues=[    max(m.getVal().tolist() )  for m in MAXMARGINALS  ]
        
        
        #MAXvalues_phred = [ PhredScore(x) for x in MAXvalues ]
        """ phred scale the error prob. so we take the 1-posterior, and phred scale it
            Note we have expoentiate the values then take the complement, then phred scale """
        phred_gq= [ round(PhredScore(x),2) for x in (1-np.exp(MAXvalues[0:totalSize])) ]
        #MAXvalues_str= [str(x) for x in MAXvalues]
        phred_str=[str(x) for x in phred_gq]
        #print "maxvalues: ",  MAXvalues[0:totalSize]
        #this is the decoding, where we get teh assignment of the variable
        # that has the max marginoal value
        MAPAssignment=MaxDecoding( MAXMARGINALS  )
        #print "MAPAssignment: ",  MAPAssignment[0:3]
        #print totalSize
        #for idx in range(totalSize):
        #    print ALPHABET[ MAPAssignment[idx] ]
        #print MAPAssignment

        #we convert from  variable assignment in the factor to genotype space
        sampleNames=pgmNetwork.returnGenotypeFactorIdxSampleNames()
      
        sampleDepths=[]
        for sample in sampleNames:
            sampleDepths.append(str(sampleDepth[sample]))
        MAPgenotypes=[ALPHABET[ i ] for i in  MAPAssignment[0:totalSize] ]
        #print MAPgenotypes

        alt=determineAltBases( MAPgenotypes, [refsequence] )
        numericalMAPgenotypes=[ numericalGenotypes(refsequence,alt,map_g) for map_g in MAPgenotypes ]
        #print numericalMAPgenotypes
        site_info="\t".join([pileup_column_chrom, str(pileup_column_pos+1), ".", refsequence,alt,qual,filter,infoString, "GT:GQ:DP"])
        #zip the genotype assignment and its log-probability together
        #genotypedata=zip(MAPgenotypes,sampleDepths,MAXvalues[0:totalSize])
        #genotypedata=zip(MAPgenotypes,sampleDepths)
        #genotypedata=zip(numericalMAPgenotypes,sampleDepths,  MAXvalues_str[0:totalSize] )
        genotypedata=zip(numericalMAPgenotypes,phred_str, sampleDepths, )
        #genotypedata=zip(numericalMAPgenotypes,sampleDepths)
        
        genotypedata=[ ":".join(list(i)) for i in  genotypedata ]
        #print genotypedata

        #print  genotypedata
        #print sampleGenotypes
        #print "\t".join( [pileup_column_chrom,  str(pileup_column_pos), refsequence]) + "\t" +"\t".join(MAPgenotypes)
        if int(pileup_column_pos) % 10000 == 0:
            sys.stdout.write("processed 10 kbp ...\n")
        #print "\t".join( [pileup_column_chrom,  str(pileup_column_pos), refsequence]) + "\t" +"\t".join(numericalMAPgenotypes)
        #print site_info
        vcfgenotypeoutput="\t".join(genotypedata)
        output=site_info+ "\t"+vcfgenotypeoutput
        vcfh.write(output+"\n")
        #if the total size (in nodes) of the network is N, the first N/2 elements
        #are genotype variables of the individuals

        #break
    
if __name__ == "__main__":
    main()
