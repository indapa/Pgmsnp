from Factor import *
from FactorOperations import *
from PedigreeFactors import *
from PgmsnpCommon import *
import itertools
#import numpy as np
import sys

"""" 
     This class is a factory for generating a genetic network
     If we consider each location in the genome independent
     we generate a new network for each position along an interval.
     A network is a collection of factors, so we return a python
     list of Factor objects. Each (read) phenotype is conditionally
     independent given its genotype, so each member of the pedigree has a
     PhenotypeGivenGenotypeFactor. If a member is a founder it is represented by
     a returnGenotypePriorFounderFactor. If a member is a non-founder it is reprsented by
     a returnGenotypeGivenParentsFactor. The first is a simple factor representing a genotype prior.
     The second is basically a Punnet square representing Mendelian inheritance probabilities.
     """

class PgmNetworkFactory(object):

    """ We construct our network based on a pedigree file, chromosome, and posiiton.
        The number of alleles is the number of possible alleles segregating at a site.
        For the full probablistic calculation it set at 4 (ACGT). In practice, assuming
        most sites are bi-allelic its value is 2."""

    def __init__(self, pedfile, chrom, position, refbase, nonfounderpriorvalues, totalAlleles=4, ploidy=2):
        #parse pedfile
        
        self.chrom=chrom
        self.pos=position
        self.refbase=refbase
        self.totalAlleles=totalAlleles
        self.pedigree=Pedfile(pedfile)
        self.pedigree.parsePedfile()
        self.pedlist=self.pedigree.getPedList()
        self.pedids=self.pedigree.returnIndivids()
        self.totalsize=self.pedigree.getTotalSize()

        #list of factors that will comprise the Genetic network
        self.totalFactors=self.pedigree.getTotalSize() * 2
        self.factorList=self.totalFactors*[None]
        self.ploidy=ploidy
        self.genotypeCardinality=len( [ "".join( list(combo)) for combo in itertools.combinations_with_replacement(['A','C', 'G', 'T'],ploidy) ] )
        self.nonfounderpriorvalues=nonfounderpriorvalues
        self.constructNetwork()



    def constructNetwork(self):
        """ make a new network  """
        totalPeople=self.pedigree.getTotalSize()
        for i in range( totalPeople ):

            j=i+1
            
            if self.pedlist[i].isFounder():
                """ if its a founder invidual, assign it GenotypePriorFounderFactor """
                self.factorList[i]=returnGenotypePriorFounderFactor( self.refbase,j )

                #previously, we assigned it a factor based on HWE and given allele frequencies.
                #I suppose we could use this for genotyping previoulsy known sites

                #self.factorList[i]=GenotypeAlleleFreqFactor(self.allelefreq,self.pedlist[i].getid(),self.pedlist[i].getid())
                #factorList(i)=genotypeGivenAlleleFreqsFactor(alleleFreqs,i);
            else:
                """ for non-founder we make GenotypeGivenParentsFactor """
                #3print self.pedlist[i].getParents(), self.pedlist[i].getid()
                #GenotypeGivenParentsFactor(2,"bart","homer","marge","""Bart | Homer, Marge """)
                #self.factorList[i]=GenotypeGivenParentsFactor(self.totalAlleles, self.pedlist[i].getid(), self.pedlist[i].getParents()[0], self.pedlist[i].getParents()[1], "child|Father,Child")
                parent1Index=self.pedids.index( self.pedlist[i].getParents()[0] )+1
                parent2Index=self.pedids.index( self.pedlist[i].getParents()[1] )+1
                child=self.pedlist[i].getid()
                
                parent1name=self.pedlist[parent1Index].getid()
                parent2name=self.pedlist[parent2Index].getid()
               
                name= " ".join ( [child+" genotype | ", "parent " + parent1name, ",",  "parent " + parent2name] )
                
                #the variable names are the index of the individuals in the list of memebers of the pedigrees
                #this helps when making the variable number for the factor Reads | Genotypes
                #self.factorList[i]=returnGenotypeGivenParentsFactor(j, parent1Index ,  parent2Index , name, self.totalAlleles)
                self.factorList[i]=returnNonFoundersFactor(j, parent1Index ,  parent2Index , self.nonfounderpriorvalues, name,self.totalAlleles)
                #returnGenotypeGivenParentsFactor(numAlleles, genotypeVarChild, genotypeVarParentOne, genotypeVarParentTwo)

            name=self.pedlist[i].getid()+" read phenotype | " + self.pedlist[i].getid() + " genotype"
            #GLFactor = Factor( [readVar, genoVar], [1,10], [], 'read_phenotype | genotype ')
            #this factor is the reads | genotype
            self.factorList[i+totalPeople]=Factor([j+totalPeople,j],[1, self.genotypeCardinality],[], name )

    def getFactorList(self):
        return self.factorList

    def printFactorList(self):
        for f in self.factorList:
            print "factor: ", f
            print 
        #return "\n".join(networkOutstrings)

    """ given a sample name, return its index position in the list of individual ids for the pgmNetwork
        If  not present, then throw a ValueError exception and exit"""
    def getSampleNamePedIndex(self,samplename):
        try:
            return self.pedids.index( samplename )
        except ValueError:
            print samplename + " not in ped file"
            sys.exit(1)


    """ return a list of the pedigree individual ids (samplenames) of the pedigree (ped file) used to construct the genetic network """
    def getSampleNames(self):
        return self.pedids


    def returnTotalSize(self):
        return self.totalsize

    def returnGenotypeFactorIdxSampleNames(self):
        """ return a list which contains the sample names  of the genotype variables in the pgm network self.factorList """
        #print range(self.totalFactors/2)
        return [ self.pedlist[i].getid()  for i in range(self.totalFactors/2) ]