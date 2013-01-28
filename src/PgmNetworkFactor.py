from Factor import *
from FactorOperations import *
from PedigreeFactors import *
import itertools
import numpy as np
"""" Still not sure how this is going to work
     This class is a factory for generating a genetic network
     If we consider each location in the genome independent
     we generate a new network for each position along an interval.
     A network is a collection of factors, so we return a python
     list of Factor objects. Each (read) phenotype is conditionally
     independent given its genotoype, so each member of the pedigree has a
     PhenotypeGivenGenotypeFactor"""

class PgmNetworkFactory(object):
    

    def __init__(self, pedfile, chrom, position, totalAlleles=4):
        #parse pedfile
        
        
        self.chrom=chrom
        self.pos=position
        self.totalAlleles=totalAlleles
        self.pedigree=Pedfile(pedfile)
        self.pedigree.parsePedfile()
        self.pedlist=self.pedigree.getPedList()
        self.pedids=self.pedigree.returnIndivids()
        

        #list of factors that will comprise the Genetic network
        self.totalFactors=self.pedigree.getTotalSize() * 2
        self.factorList=self.totalFactors*[None]

        
    def constructNetwork(self):
        totalPeople=self.pedigree.getTotalSize()
        for i in range( totalPeople ):
            
            
            if self.pedlist[i].isFounder():
                self.factorList[i]=Factor()
                #print self.pedlist[i].getid()

                self.factorList[i]=Factor( [self.pedlist[i].getid()], [], [], self.pedlist[i].getid())
                #self.factorList[i]=GenotypeAlleleFreqFactor(self.allelefreq,self.pedlist[i].getid(),self.pedlist[i].getid())
                #factorList(i)=genotypeGivenAlleleFreqsFactor(alleleFreqs,i);
            else:
                #3print self.pedlist[i].getParents(), self.pedlist[i].getid()
                #GenotypeGivenParentsFactor(2,"bart","homer","marge","""Bart | Homer, Marge """)
                #self.factorList[i]=GenotypeGivenParentsFactor(self.totalAlleles, self.pedlist[i].getid(), self.pedlist[i].getParents()[0], self.pedlist[i].getParents()[1], "child|Father,Child")
                parent1Index=self.pedids.index( self.pedlist[i].getParents()[0] )
                parent2Index=self.pedids.index( self.pedlist[i].getParents()[1] )
                child=self.pedlist[i].getid()
                parent1name=self.pedlist[parent1Index].getid()
                parent2name=self.pedlist[parent2Index].getid()
                name=child+" genotype |"+parent1name+","+parent2name
                self.factorList[i]=returnGenotypeGivenParentsFactor(self.totalAlleles, i, parent1Index ,  parent2Index , name)
                #returnGenotypeGivenParentsFactor(numAlleles, genotypeVarChild, genotypeVarParentOne, genotypeVarParentTwo)
            name=self.pedlist[i].getid()+" read phenotype | " + self.pedlist[i].getid() + " genotype"
            self.factorList[i+totalPeople]=Factor(i+totalPeople,i, name )

    def getFactorList(self):
        return self.factorList


