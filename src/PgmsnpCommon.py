#!/usr/bin/env python
from Genotype import *
import numpy as np
import itertools
import sys
from common import *
from Factor import *
from PGMcommon import *
from FactorOperations import *

def calculateGLL(pileup):
    """ calculate genotype log-likelhood(s) for all 10 possible genotypes
        given a list of namedtuples representing a pileup column """
    genotypes=[] #list of Genotype objects
    ploidys=[2]
    #enumerate possible genotype combinations for the value of ploidy
    #with the bases A,C,G,T
    for p in ploidys:
        l=[ "".join( list(combo)) for combo in itertools.combinations_with_replacement(['A','C','G','T'],p) ]
        for g in l:
            #genotype="".join( list(g) )
            genotypes.append( Genotype( g, p) )

    h=len(pileup)
    #print "total pileup: ", h
    #the likelihood matrix is a 10 column matrix
    #the number of rows is equal to the pileup count i.e. coverage at position
    #we pre-allocate the likelihood matrix
    likelihood_matrix=np.zeros( (h,10) )
    #now we iterate through pileup, taking the ith called base and quality
    for i in range(h):
        #print pileup[i]
        (calledBase, phred)=(pileup[i].basecall, pileup[i].bq)
        #for the jth possible genotype
        for j in range( len( genotypes)):
            #get the error probability
            errorprob=ErrorProb(phred)
            #calculate the base likelihood (Pr(B=b|G=g)
            #fill in likelihood for the ith pileuped base for the jth possible genotype
            likelihood_matrix[i,j]=genotypes[j].calculateBaseLikelihood(calledBase,errorprob)


    #now this is the log likelihood 
    #print np.log(likelihood_matrix)
    #the overall (log) likelihood are each of the 10 column sums
    #genotypelikelihoods=np.sum( np.log(likelihood_matrix), axis=0)
    return np.sum( np.log(likelihood_matrix), axis=0)


def indexToGenotype( index, ploidy=2,alleles='ACGT' ):
    """ return genotype at a given index position after
    enumerating all possible genotypes given string of alleles and
    assigning to a list. By default the list contains all possible 10 genotypes"""


    genotypes= [ "".join(list(genotype))  for genotype in itertools.combinations_with_replacement(alleles, ploidy) ]

    try:
       return genotypes[index]
    except IndexError:
        print "Index out of bounds, not a valid index for list of genotypes"

def genotypeToIndex( geno, ploidy=2,alleles='ACGT'):
    """ given a genotype return its index in the enumerated list of possible genotypes
    given ploidy and alleles """


    genotypes= [ "".join(list(genotype))  for genotype in itertools.combinations_with_replacement(alleles, ploidy) ]
    

    try:
        return  genotypes.index(geno)
    except ValueError:
        print "genotype not in list of genotypes."


def returnGenotypePriorFounderFactor( refbase,altbases,genotypePrior, theta=0.001,ploidy=2):
    """ Not sure this is right, but its simple enough
        This function returns a factor representing genotype priors, passing in the
        reference base, and list of alt alelles in altbase. genotypePrior is the name of the variable
        theta is heterozygosity rate set to .001 by default and ploidy is set to 2
        prior(ref homoz=1-3(theta/2), het=theta, alt homoz=theta/2  """

    alleleList=([refbase]+altbases)
    print ''.join(alleleList)
    numAlleles=len(alleleList)
    f1= Factor( [genotypePrior ], [ ], [ ], 'genotypePrior')
    (allelesToGenotypes, genotypesToAlleles)=generateAlleleGenotypeMappers(numAlleles)
    (ngenos,ploidy)=np.shape(genotypesToAlleles)
    f1.setCard([ ngenos] )
    values=np.zeros( (np.prod(f1.getCard()))).tolist()
    #print values
    l=[ "".join( list(combo)) for combo in itertools.combinations_with_replacement(alleleList,ploidy) ]
    #l=[ "".join( list(combo)) for combo in itertools.product(alleleList,repeat=ploidy) ]
    print l
    for i in range(ngenos):
        genotype=indexToGenotype(i, ''.join(alleleList) )
        (a1,a2)=list(genotype)
        #print a1,a2
        if a1 == a2 and refbase not in genotype:
            print genotype, 'non-ref homoz'
            values[i]=(theta / 2.)
        elif a1==a2==refbase:
            print genotype, 'homoz reference'
            values[i]= 1 - (3*(theta/2.))
        elif a1!=a2 and refbase in genotype:
            print genotype, 'heterzygote'
            values[i]=theta
        else:
            print genotype, 'tri-alleleic het'
            values[i]=np.power( [ theta/2 ], 3).tolist()[0]
    #print values
    f1.setVal(values)
    return f1

def returnGenotypeGivenParentsFactor( numAlleles, genotypeVarChild, genotypeVarParentOne, genotypeVarParentTwo  ):
    """ return a Factor object that represents pr( offspring_genotype | genotype_mother, genotype_father )
        basically this is a Punnet square """
    f1= Factor( [genotypeVarChild, genotypeVarParentOne, genotypeVarParentTwo ], [ ], [ ])
    (allelesToGenotypes, genotypesToAlleles)=generateAlleleGenotypeMappers(numAlleles)
    (ngenos,ploidy)=np.shape(genotypesToAlleles)
    f1.setCard([ ngenos,ngenos,ngenos ] )
    #set the values to zero initially
    values=np.zeros( (np.prod(f1.getCard()))).tolist()
    assignments=IndexToAssignment( np.arange(np.prod(f1.getCard())), f1.getCard() )-1
    for z in range( np.prod(f1.getCard() ) ):
        curr_assign= assignments[z]
        childAssignment=int(curr_assign[0])

        parent1gametes= genotypesToAlleles[curr_assign[1],:]
        parent2gametes= genotypesToAlleles[curr_assign[2],:]
        #print 'parental gametes: ', parent1gametes, parent2gametes
        #print 'child assignment: ', childAssignment
        #list of tuples containing list of zygote(genotype) tuples
        zygote_list=list(itertools.product(parent1gametes,parent2gametes))
        punnet_freq=[  allelesToGenotypes[zygote[0],zygote[1]] for zygote in zygote_list ]
        histc={}
        hist=[]
        for g in range( ngenos):
            histc[g]=0.
        for x in punnet_freq:
            histc[x]+=1.
            #print histc.values()
        for g in range (ngenos):
            hist.append ( histc[g] )
            #print punnet_freq
        hist=(np.array ( hist)) /4
        values[z]=hist[childAssignment]

        f1.setVal( values )

    return f1
