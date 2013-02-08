import numpy as np
class Genotype(object):
    """ represents a genotype can be of arbritary ploidy from 1 ... N
        Genotype is unphased e.g. AC would be the same CA"""

    def __init__(self, genotype='NN', ploidy=2):
        """ default ploidy is 2, default genotype is NN """

        self.genotype=genotype
        self.ploidy=ploidy

    def __str__(self):
        return "\t".join( [self.genotype, str(self.ploidy)])

#def __cmp__(self,other):
#    """ need to implement
#        comparison based on genotype"""


    def assertGenotype(self):
        assert (len(self.genotype) == self.ploidy), "ploidy of %d incompatiable with genotype %s" % (self.ploidy, self.genotype)

    def setPloidy(self, ploidy):
        """ set ploidy of gentoype """
        self.ploidy=ploidy

    def getPloidy(self):
        """ get ploidy of genotype """
        return self.ploidy

    def isHomoz(self):
        return all(x == items[0] for x in items)

    """ given a specific base and associated base error probbablility,
    calculate Pr(base | self.genotype ). Note, the base error probablity is not in Phred space
    for more see http://en.wikipedia.org/wiki/Phred_quality_score#Reliability
    """
    def calculateBaseLikelihood ( self, specificBase, baseErrorProb):
        psGenotype=0.
    
        self.assertGenotype() #check if the ploidy is compatabile with self.genotype
        for i in range(self.ploidy):
            if self.genotype[i] != '-':
                if specificBase == self.genotype[i]:
                    psGenotype += ( 1 - baseErrorProb )
                else:
                    psGenotype += baseErrorProb/3


        return psGenotype / self.ploidy


    """ The genotype likelihood (GL) is defined as Pr( read | G )
            Let N be the total number of mapped/pileup bases.
            N=N_a + N_c + N_g + N_t where each term is count of each
            pileuped base.
            For example if G(enotype) is AA. The GL is defined as
            Pr(read|AA) = prod_j=1^N_a = (1-e_j) * prod_k=1^N-N_a ( e_k/3 )
            where e_i or e_j is the base error probability. We divide by three
            assuming a base and equal prob. of being miscalled as any of the 3 alternative
            bases. """

    def calculateCorrectGenotypeLikelihood ( self, read, errorprob ):

        psGenotype=0.

        self.assertGenotype() #check if the ploidy is compatabile with self.genotype

        if all(x == self.genotype[0] for x in self.genotype):
            #print "homoz ", self.genotype

            
            if read == self.genotype[0]:
                psGenotype =  1 - errorprob 
            else:
                psGenotype =  errorprob/3 


        else:
            #print "het ", self.genotype

            if read in self.genotype:
                psGenotype =  .5 * (1 - (2.*errorprob)/3  )  
            else:
                psGenotype = errorprob/3 

        return psGenotype
            

