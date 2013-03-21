import sys
class Bedfile(object):
    """ represent a bedfile """

    def __init__(self, bedfilename):
        self.name=bedfilename
        self.fh=None

    def open(self):
        try:
            self.fh=open(self.name, 'r')

        except IOError:
            sys.stderr.write("Error: cannot open bedfile for reading.\n")

    def yield_bedcoordinate(self):
        """ yield a tuple of (chr, start,end) from bed file """
        for line  in self.fh:
            if '@' in line: continue
            fields=line.strip().split("\t")
            (chr, start, end) = fields[0:3]
            yield(chr, int(start), int(end) )
