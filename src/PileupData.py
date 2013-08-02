import itertools
import collections
class PileupData(object):
    """ represents a pileup column of reads for a chromosome and posiiton
        please note, start/stop is zero based, half open interval consistent
        with bed file"""


    def __init__(self, chrom, start, end, pos, pileuplist=[]):
        self.chrom=chrom
        self.start=start
        self.end=end
        self.pileupcolumn_pos=pos
        """ this will be a list of tuples
            with the tuples consisting of (sample,readgroup,alignment name, base, base quality"""
       
        Pileup = collections.namedtuple('Pileup', ['sample', 'RG', 'alignmentname', 'basecall', 'bq'])
        #self.pileupList is now a list of Pileup namedtuple(s)
        self.pileupList= map(Pileup._make, pileuplist)

   
    def __str__(self):
        return "\t".join( [self.chrom, str(self.start), str(self.end), str(self.pileupcolumn_pos), "total pileupd reads: " + str(len(self.pileupList))])

    def getBedString(self):
        return "\t".join([self.chrom,str(self.start), str(self.end) ])

    def getChrom(self):
        return self.chrom
    def setChrom(self,chrom):
        self.chrom=chrom

    def getStart(self):
        return self.start
    def setStart(self,start):
        self.start=start

    def getEnd(self):
        return self.end
    def setEnd(self,end):
        self.end=end

    def getCoverage(self):
        """ the coverage is the total number of tuples  in the pileupList list """
        return len(self.pileupList) 
    
    def getPileupColumnPos(self):
        """ return the pileup column position """
        return self.pileupcolumn_pos

    def yieldSamplePileupData(self):
        """ we can yield a list of pileup data by grouping info by sample
            see http://docs.python.org/release/3.2.3/library/itertools.html#itertools.groupby
            and http://stackoverflow.com/a/7286/1735942 on how to use itertools.grouper method
            for pileupdata in PileupOBject.yieldSamplePileupData:"""
        
        # list needs to be sorted: http://stackoverflow.com/a/12088018/1735942
        self.pileupList = sorted(self.pileupList, key=lambda x: x.sample)
        #http://stackoverflow.com/a/783/1735942
        #groupby needs the dat to be sorted becaue groupby method actually just iterates through a list and whenever the key changes it creates a new group.
        for key, group in itertools.groupby(self.pileupList, lambda x: x.sample ):
            yield ( key,list(group) )

    def getSamplePileupData(self, samplename):
        """ return list of namedtuple of pileup data for the particular samplename   """
        
        #return [ tple for tple in self.pileupList if tple[0]== samplename ]
        return [ tple for tple in self.pileupList if tple.sample== samplename ]

    def getPileupDepth(self):
        """ return the total number of reads piled up """
        return len(self.pileupList)

    def getSetObservedBases(self):
        """ return a set object of observed bases for the entire pileup of bases """
        return set( [ namedtpl.basecall for namedtpl in self.pileupList ] )

    

