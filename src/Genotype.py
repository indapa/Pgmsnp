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

    """ given a specific base and associated base error probbablility,
    calculate Pr(self.genotype | base). Note, the base error probablity is not in Phred space
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