#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/14/16 4:59 PM
# @Author  : Runsheng
# @Email   : Runsheng.lee@gmail.com
# @File    : test_primer_check.py


import unittest

# self import
from primervcf.map2vcf import *
from primerdiffer.utils import dic2dic, fasta2dic

class TestPrimerCheck(unittest.TestCase):

    def setUp(self):
        self.cb4_genome = "/t1/ref_BN/cb5.fa"
        self.cb4=dic2dic(fasta2dic(self.cb4_genome))

    def test_map2vcf(self):
        # use the file in fastq-dump --split-spot SRR1793006
        pass


    def tearDown(self):
        self=None

if __name__ == '__main__':
    unittest.main()