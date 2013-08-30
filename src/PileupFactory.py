from PileupData import PileupData
#from BedFile import Bedfile
import pdb
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




    def yieldPileupData(self, minmapq=0, minbq=0):
        """ iterate thru the bedobj bed intervals
            processing the pileupcolumns
            for each of the reads in the pileupcolumn get the RG, SM, basecall, and base quality"""
        #here we just iterate thru the pileup columns collecting the observed basecalls and qualities in list of tuples

        for (chrom, start, end) in self.bedobj.yield_bedcoordinate():
            

            """ the pileup engine gets all reads that overlap the bed interval
                if the pileupcolumn.pos is out of range of the bed interval, skip it"""

            pileup_iter=self.pybam.pileup( chrom,start,end,stepper = "all" )
            for pileupcolumn in pileup_iter:
                
                if pileupcolumn.pos+1 > end: continue
                if pileupcolumn.pos < start: continue

                observed_data=[] #a list of (basecall, basequality) tuples
                for pileupread in pileupcolumn.pileups:
                    #pdb.set_trace()
                    if pileupread.alignment.mapq < minmapq: continue
                    
                    """ skip deletions and insertions for now ... """
                    if pileupread.is_del or pileupread.indel !=0: 
                        #pdb.set_trace()
                        continue
                    
                    """ skip N bases """
                    if pileupread.alignment.seq[pileupread.qpos] == 'N':continue
                    
                    bq=ord ( pileupread.alignment.qual[ pileupread.qpos ] )  - 33
                    if bq < minbq: continue
                   
                    #if bq == 0: continue
                    #print pileupread.alignment
                    #tid=pileupread.alignment.tid
                   
                    readgroup=dict( pileupread.alignment.tags )['RG']
                    
                    sample=self.readgroupdict[readgroup]
                    #print 'mapq: ', pileupread.alignment.mapq, " bq: ", bq, pileupread.alignment.qname, sample
                    #ascii conver the basequality
                    
                    
                    #print samfile.getrname(tid),pileupcolumn.pos, observed_data
                    #observed_data is a list of tuples
                    #the tuple is [samplename, readgroupID, alignment_name,basecall,phred_scale_basequality
                    observed_data.append( (sample, readgroup, pileupread.alignment.qname, pileupread.alignment.seq[pileupread.qpos], bq, pileupread.alignment.cigarstring) )
                yield ( PileupData( chrom, start, end, pileupcolumn.pos, observed_data ) )

                    


