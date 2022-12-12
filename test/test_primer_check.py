#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/14/16 4:59 PM
# @Author  : Runsheng
# @Email   : Runsheng.lee@gmail.com
# @File    : test_primer_check.py


import unittest,os

# self import
from primervcf.primer_indel import *
from primerdiffer.utils import dic2dic, fasta2dic

class TestPrimerCheck(unittest.TestCase):

    def setUp(self):
        self.sp9_genome = "/t1/ref_BN/cn3_new.fa"
        self.cb4_genome = "/t1/ref_BN/cb5.fa"
        self.cb4=dic2dic(fasta2dic(self.cb4_genome))

    def test_parse_cigar(self):
        print(parse_cigar(
            "75S2M1I76M2I8M1D35M1I2M1I38M1D21M"))

    def test_get_del(self):
        wkdir="/t1/ref_BN/vcf/SRR1793006"
        os.chdir(wkdir)
        vcf_file="QR25.vcf"

        QR25=get_deletion_region(vcf_file, key="del",len_cutoff=10)
        print(len(QR25))
        print(QR25[1:10])
        bedl2file(QR25, "QR25_del.bed")

    def test_design_primer_one(self):
        wkdir="/t1/ref_BN/vcf/SRR1793006"
        db = "/t1/ref_BN/cb5.fa"
        os.chdir(wkdir)
        vcf_file="QR25_del.bed"
        QR25=file2bedl(vcf_file)
        p=del_primer_check(QR25[2], db=db, ref_dict=self.cb4, debugmod=True)
        print(p)

    def test_walk_deletion(self):
        wkdir="/t1/ref_BN/vcf/SRR1793006"
        db = "/t1/ref_BN/cb5.fa"
        os.chdir(wkdir)
        vcf_file="QR25_delhead.bed"
        QR25=file2bedl(vcf_file)
        primer_d=flow_walk_deletion(ref_dict=self.cb4, db=db, deletion_bedlist=QR25, prefix="qr25")
        print(primer_d)


    def tearDown(self):
        self=None

if __name__ == '__main__':
    unittest.main()