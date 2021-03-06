#!/usr/bin/env python
from Genotype import *
import numpy as np
import itertools
from common import *
from Factor import *
from PGMcommon import *
from FactorOperations import *
from collections import defaultdict
import pdb

def calculateGLL(pileup,ploidy=2):
    """ calculate genotype log-likelhood(s) for all 10 possible genotypes
        given a list of namedtuples representing a pileup column """
    genotypes=[] #list of Genotype objects
    
    """ enumerate possible genotype combinations for the value of ploidy with the bases A,C,G,T """
   
    l=[ "".join( list(combo)) for combo in itertools.combinations_with_replacement(['A','C','G','T'],ploidy) ]
    #print l
    genotypes=[ Genotype( g, ploidy)  for g in l  ]

    #for g in l:
        
    #    genotypes.append( Genotype( g, ploidy) ) #make a Genotype object, add to genotypes list

    total_genotypes=len(genotypes)
    h=len(pileup)
    #print "total pileup: ", h
    """ the likelihood matrix is a 10 column matrix, where columns are genotypes
       the number of rows is equal to the pileup count i.e. coverage at position
       we pre-allocate the likelihood matrix"""
    likelihood_matrix=np.zeros( (h,total_genotypes) )
    #now we iterate through pileup, taking the ith called base and quality """
    for i in range(h):
        #print pileup[i]
        (calledBase, phred)=(pileup[i].basecall, pileup[i].bq)
        #print calledBase
        #for the jth possible genotype
        for j in range( len( genotypes)):
            
            #get the error probability
            errorprob=ErrorProb(phred)
            #calculate the base likelihood (Pr(B=b|G=g)
            #fill in likelihood for the ith pileuped base for the jth possible genotype
            #likelihood_matrix[i,j]=genotypes[j].calculateBaseLikelihood(calledBase,errorprob)
            likelihood_matrix[i,j]=genotypes[j].calculateCorrectGenotypeLikelihood(calledBase,errorprob)
    """ debug print statments """
    #for g in genotypes: print g
    #print np.shape( likelihood_matrix )
    #print np.log(likelihood_matrix)

    """now this is the log likelihood """ 
    #print np.log(likelihood_matrix)
    """the overall (log) likelihood are each of the 10 column sums
       we return the the genotype (log) likelihood for each possible genotype """
    
    
    return np.sum( np.log(likelihood_matrix), axis=0)


def calculateNoObservationsGL(ploidy=2):
    """ the numbers returned are in probability space, not log-space! """
    
    l=[ "".join( list(combo)) for combo in itertools.combinations_with_replacement(['A','C','G','T'],ploidy) ]
    genotypes=[ Genotype( g, ploidy)  for g in l ]
    total_genotypes=len(genotypes)
    """ the prob(no data | genotype ) is 1
        We keep the numbers in probablity space"""
    likelihood_matrix=np.ones( (1,total_genotypes))
    return np.sum( likelihood_matrix, axis=0)


def indexToGenotype( index,alleles='ACGT',ploidy=2 ):
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


def returnGenotypePriorFounderFactor( refbase, factorVar, theta=0.001,ploidy=2 ):
    """ Not sure this is right, but its simple enough
        This function returns a factor representing genotype priors, passing in the
        reference base, and list of alt alelles in altbase. genotypePrior is the name of the variable
        theta is heterozygosity rate set to .001 by default and ploidy is set to 2
        prior(ref homoz=1-3(theta/2), het=theta, alt homoz=theta/2  """

    
    numAlleles=len( ['A','C','G','T'] )
    f1= Factor( [factorVar ], [ ], [ ], 'genotypePrior')
    (allelesToGenotypes, genotypesToAlleles)=generateAlleleGenotypeMappers(numAlleles)
    (ngenos,ploidy)=np.shape(genotypesToAlleles)
    #print ngenos
    f1.setCard([ ngenos] )
    values=np.zeros( (np.prod(f1.getCard()))).tolist()
    #print values
    # l is the exhaustive set possible genotypes for a given ploidy
    #l=[ "".join( list(combo)) for combo in itertools.combinations_with_replacement(['A','C','G','T'],ploidy) ]
    #print l
    
    for i in range(ngenos):
        
        genotype=indexToGenotype(i, ''.join( ['A','C','G','T'] ) )
        (a1,a2)=list(genotype)
        #print a1,a2
        if a1 == a2 and refbase not in genotype:
            #print genotype, 'non-ref homoz'
            values[i]=(theta / 2.)
        elif a1==a2==refbase:
            #print genotype, 'homoz reference'
            values[i]= 1 - (3*(theta/2.))
        elif a1!=a2 and refbase in genotype:
            #print genotype, 'heterzygote'
            values[i]=theta
        else:
            #print genotype, 'tri-alleleic het'
            values[i]=np.power( [ theta/2 ], 3).tolist()[0]
    #print values
    f1.setVal(values)
    return f1


def returnNonFoundersFactor(  genotypeVarChild, genotypeVarParentOne, genotypeVarParentTwo, values, factorName="child|parent 1, parent2",numAlleles=4  ):
    """ return a Factor object that represents pr( offspring_genotype | genotype_mother, genotype_father )
        values are the transition probalities of pr(offspring_genotype|mother,father)
        These don't change, so we calculate them once and the pass them in as a parameter.
        The only thing you are doing is setting the variable names and cardinality (based on the number of alleles)
        Note, when you calculate the transition probablities in values with returnPunnetValues,
        make sure the numAlleles is the same. Otherwise there will be a dimenionality mismatch!"""
    f1= Factor( [genotypeVarChild, genotypeVarParentOne, genotypeVarParentTwo ], [ ], values, factorName )
    (allelesToGenotypes, genotypesToAlleles)=generateAlleleGenotypeMappers(numAlleles)
    (ngenos,ploidy)=np.shape(genotypesToAlleles)
    f1.setCard([ ngenos,ngenos,ngenos ] )
    #set the values to zero initially
    
    

    return f1


def returnGenotypeGivenParentsFactor(  genotypeVarChild, genotypeVarParentOne, genotypeVarParentTwo, factorName="child|parent 1, parent2", numAlleles=4  ):
    """ return a Factor object that represents pr( offspring_genotype | genotype_mother, genotype_father )
        basically this is a Punnet square """
    f1= Factor( [genotypeVarChild, genotypeVarParentOne, genotypeVarParentTwo ], [ ], [ ], factorName )
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
        histc = defaultdict(int)
        hist=[]
        
        for x in punnet_freq:
            histc[x]+=1.
            #print histc.values()

        hist=[ histc[g] for g in range(ngenos) ]
        #for g in range (ngenos):
        #    hist.append ( histc[g] )
            #print punnet_freq
        hist=(np.array ( hist)) /4
        values[z]=hist[childAssignment]

        f1.setVal( values )

    return f1


def returnPunnetValues( numAlleles=4 ):
    """ return a list of transition probabilities for an offpsring genotype, conditional
        on parental genotypes: Pr(offspring| parent1, parent2)
        Basically this returns the probabilities from  Punnet square.
        We pass in the number of alleles at a site. If we consider a bi-allelic site,
        then we return a 27 transition probabilities, since there are 3x3x3 possible genotypes
        given 2 alleles.

        Complete enumeration be passing in 4 alleles, for a possible of 10 genotypes for each
        individual. This would be a 10x10x10 size of transition probabilities.

        Father Mother offspring_homoz offspring_het offspring_het
        dd       dd          1              0             0
        dd       Dd         .5             .5             0
        .
        .
        .
        DD       DD          0              0             0"""

    (allelesToGenotypes, genotypesToAlleles)=generateAlleleGenotypeMappers(numAlleles)
    (ngenos,ploidy)=np.shape(genotypesToAlleles)

    cardinality=np.array( [ngenos,ngenos,ngenos] )

    #set the values to zero initially
    
    values=np.zeros( np.prod(cardinality)).tolist()
    assignments=IndexToAssignment(np.arange(np.prod(cardinality)),cardinality)-1

    for z in range( np.prod(cardinality) ):
        curr_assign= assignments[z]

        childAssignment=int(curr_assign[0])

        parent1gametes= genotypesToAlleles[curr_assign[1],:]
        parent2gametes= genotypesToAlleles[curr_assign[2],:]

        #list of tuples containing list of zygote(genotype) tuples
        zygote_list=list(itertools.product(parent1gametes,parent2gametes))

        punnet_freq=[  allelesToGenotypes[zygote[0],zygote[1]] for zygote in zygote_list ]

        histc = defaultdict(int)

        hist=[]

        for x in punnet_freq:
            histc[x]+=1.
            #print histc.values()
        hist=[ histc[g] for g in range(ngenos) ]


        hist=(np.array ( hist)) /4
        values[z]=hist[childAssignment]
        if values[z] == 0:
            values[z] = .00001

        #print

    return values



def generateAllReferenceGenotypesAssignment ( g_idx, totalGenotypes):
    """ return the assignment that corresponds to all individuals being homozygous reference
    g_idx is the index of the reference genotype in the samplespace of genotypes returned by
    calling genotypeToIndex.
    totalGenotypes is the number genotypes which is half the number of total variables in the pgm network """
    total=totalGenotypes * 2
    assignment= [g_idx for i in range(0, totalGenotypes) ]
    for i in range(totalGenotypes, total):
        assignment.append(0)
    return np.array (assignment)






