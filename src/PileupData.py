from itertools import groupby
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
            with the tuples consisting of (readgroup,sample,alignment name, base, base quality"""
        self.pileupList=pileuplist

   
    def __str__(self):
        return "\t".join( [self.chrom, str(self.start), str(self.end), str(self.pileupcolumn_pos)])

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
        """ the coverage is the total number of records  in the piluepcolumn_pos list """
        return len( self.pileupcolumn_pos)

    def yieldSamplePileupData(self):
        """ we can yield a list of pileup data by grouping info by sample
            see http://docs.python.org/release/3.2.3/library/itertools.html#itertools.groupby
            and http://stackoverflow.com/a/7286/1735942 on how to use itertools.grouper method

            for pileupdata in PileupOBject.yieldSamplePileupData:"""
        for key, group in groupby(self.pileupList, lambda x: x[1]):
            yield list(group)

    def getSamplePileupData(self, samplename):
        """ return the pileupdata for the particular samplename   """
        return [ tple for tple in self.pileupList if tple[0]== samplename ]

    

