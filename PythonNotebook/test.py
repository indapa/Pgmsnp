import unittest
import pysam
import sys
from BedFile import Bedfile
from PileupData import PileupData
from PileupFactory import PileupFactory
from PgmsnpCommon import *
from common import *
import pdb

DATADIR="/Users/amit/data/MySimulations/Simulation1"
PEFILE=DATADIR+"/fam1.ped"
BEDFILE="/Users/amit/software/Pgmsnp/PythonNotebook/test.sim.bed"
BAMFILE=DATADIR+"/fam1.bed"

class PgmsnpTest(unittest.TestCase):
    
    def setUp(self):
        self.pybamfile=pysam.Samfile( BAMFILE, 'rb')
        
    def testPileup(self):
        bedobj=Bedfile(BEDFILE)
        Pfactory=PileupFactory(pybamfile,bedobj)
