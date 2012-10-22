from PileupData import PileupData
#from BedFile import Bedfile

class PileupFactory(object):
    """ generates PileupData objects given a Pysam and bed file object """

    def __init__(self, pysamobj,bedfileobj):
        self.pybam=pysamobj
        self.bedobj=bedfileobj
        self.bedobj.open()
        self.readgroupdict={}
        self.initializeRGdict()

    def initializeRGdict(self):
        rgdictlist=self.pybam.header['RG']
        for dictionary in rgdictlist:
            #print dictionary['ID'], dictionary['SM']
            self.readgroupdict[ dictionary['ID'] ]= dictionary['SM']




    def yieldPileupData(self):
        """ iterate thru the bedobj bed intervals
            processing the pileupcolumns
            for each of the reads in the pileupcolumn get the RG, SM, basecall, and base quality"""
        #here we just iterate thru the pileup columns collecting the observed basecalls and qualities in list of tuples

        for (chrom, start, end) in self.bedobj.yield_bedcoordinate():
            
            pileup_iter=self.pybam.pileup( chrom,start,end,stepper = "all" )
            for pileupcolumn in pileup_iter:

                observed_data=[] #a list of (basecall, basequality) tuples
                for pileupread in pileupcolumn.pileups:
                    #print pileupread.alignment.qname
                    tid=pileupread.alignment.tid
                    readgroup=dict( pileupread.alignment.tags )['RG']
                    sample=self.readgroupdict[readgroup]
                    #ascii conver the basequality
                    bq=ord ( pileupread.alignment.qual[ pileupread.qpos ] )  - 33
                   
                    #print samfile.getrname(tid),pileupcolumn.pos, observed_data
                    observed_data.append( (sample, readgroup, pileupread.alignment.qname,pileupread.alignment.seq[pileupread.qpos], bq) )
                yield ( PileupData( chrom, start, end, pileupcolumn.pos, observed_data ) )

                    


