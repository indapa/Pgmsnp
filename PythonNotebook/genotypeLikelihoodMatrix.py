# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

# <codecell>

""" calculate the genotype likelihood for all possible 10 genotypes
    given the list of tuples [ (basecall, phred-score) ] """

from Genotype import *
from common import *
import numpy as np
from itertools import *
genotypes=[] #list of Genotype objects
ploidys=[2]
#enumerate possible genotype combinations for the value of ploidy
#with the bases A,C,G,T
for p in ploidys:
    l=[ combo for combo in combinations_with_replacement(['A','C','G','T'],p) ]
    for g in l:
        genotype="".join( list(g) )
        genotypes.append( Genotype( genotype, p) )

pileup=[('A', 22), ('A', 4), ('T', 29)]
h=len(pileup)
#the likelihood matrix is a 10 column matrix
#the number of rows is equal to the pileup count i.e. coverage at position
#we pre-allocate the likelihood matrix
likelihood_matrix=np.zeros( (h,10) )

#now we iterate through pileup, taking the ith called base and quality
for i in range(h):
    (calledBase, phred)=pileup[i]
    #for the jth possible genotype
    for j in range( len( genotypes)):
        #get the error probability
        errorprob=ErrorProb( phred)
        #calculate the base likelihood (Pr(B=b|G=g)
        #fill in likelihood for the ith pileuped base for the jth possible genotype
        likelihood_matrix[i,j]=genotypes[j].calculateBaseLikelihood(calledBase,errorprob)


#now this is the log likelihood 
#print np.log(likelihood_matrix)
#the overall (log) likelihood are each of the 10 column sums
genotypelikelihoods=np.sum( np.log(likelihood_matrix), axis=0)

# <codecell>


