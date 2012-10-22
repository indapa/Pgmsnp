# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

# <codecell>

import pysam
import sys
from BedFile import Bedfile
from PileupData import PileupData
from PileupFactory import PileupFactory

""" given a BAM file and Bedfile, iterate through pileup columns defined by the bed interval """

bamfilename='/Users/amit/BCResearch/OlivierData/Project_TOPSWga/AlignmentGazing/wga.SNP.novel.12:82756:82957.bam'
pybamfile = pysam.Samfile(bamfilename, "rb" )
bedobj=Bedfile('/Users/amit/software/Pgmsnp/src/test.bed')
Pfactory=PileupFactory(pybamfile,bedobj)

#tell the factory object to make PileupData objects
for pileup_obj in Pfactory.yieldPileupData():
    print pileup_obj
    # the PileupData object contains pileup for all samples in the bam
    # we can iterate through pileup column by each sample's pileup
    for pileup_sample_data in pileup_obj.yieldSamplePileupData():
        print pileup_sample_data

# <codecell>


