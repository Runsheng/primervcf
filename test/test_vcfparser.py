#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/14/16 4:59 PM
# @Author  : Runsheng
# @Email   : Runsheng.lee@gmail.com
# @File    : test_primer_check.py


import unittest,os

# self import
from primerdiffer.utils import dic2dic, fasta2dic

from primervcf.vcfparser import parse_cigar, get_deletion_region, bedl2file, file2bedl


class TestPrimerCheck(unittest.TestCase):

    def setUp(self):
        self.sp9_genome = "/t1/ref_BN/cn3_new.fa"
        self.cb4_genome = "/t1/ref_BN/cb5.fa"
        self.cb4=dic2dic(fasta2dic(self.cb4_genome))

    def test_parse_cigar(self):
        print(parse_cigar(
            "75S2M1I76M2I8M1D35M1I2M1I38M1D21M"))

    def test_get_del(self):
        # self NGS based indels
        wkdir="/t1/ref_BN/vcf/SRR1793006"
        os.chdir(wkdir)
        vcf_file="QR25.vcf"
        QR25=get_deletion_region(vcf_file, key="del",len_cutoff=10)
        print(len(QR25))
        print(QR25[1:10])
        bedl2file(QR25, "QR25_del.bed")

    def test_get_del_svcaller(self):
        # caller from SV do not use indel, but use INS and DEL
        wkdir=("./SV")
        os.chdir(wkdir)
        vcf_file = "sv.vcf"
        sv= get_deletion_region(vcf_file, key="del", len_cutoff=10)
        print(len(sv))
        bedl2file(sv, "sv_del.bed")



    def tearDown(self):
        self=None

if __name__ == '__main__':
    unittest.main()