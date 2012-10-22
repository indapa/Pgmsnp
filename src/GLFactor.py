import sys
import numpy as np
from Factor import *
from PGMcommon import *
from common import *
from FactorOperations import *
from PileupData import *
from Genotype import *
from itertools import *


class GenotypeLikelihoodFactor(object):
    """ Represent the conditional probablity distribution of Pr( basecall | genotype )
    The values would represent the genotype data likelihood of seeing a set of
    basecalls in a pileup column of a BAM/SAM file assuming a genotype"""

    def __init__(self, pileup_data, pileup_var, samplename):
        
        self.pileup=pileup_data
        self.name=name
        self.genotypeList=[]
        for ploidy in [2]:
            l=[ combo for combo in combinations_with_replacement(['A','C','G','T'],p) ]
            for g in l:
                genotype="".join( list(g) )
                self.genotypeList.append( Genotype( genotype, ploidy) )


        GL=Factor ( [ pileup_var ],[ 10 ], [], samplename )

        depth=length(pileup_data)
        likelihood_matrix=np.zeros(  depth, 10 )

        for i in range(depth):
            (sample, readgroup, aligned_read, basecall, bq)=self.pileup[i]
            for j in range( len (self.genotypeList)):
                likelihood_matrix[i,j]=self.genotypeList[j].calculateBaseLikelihood(basecall, ErrorProb( bq) )
        genotypelikelihoods=np.sum( np.log(likelihood_matrix), axis=0)
        GL.setVal(genotypelikelihoods.tolist())

        """ https://github.com/indapa/Pgmsnp/blob/master/PythonNotebook/genotypeLikelihoodMatrix.py  """

        

