#!/usr/bin/env python
#-*- coding: utf-8 -*-
# @Time    : 12/12/2022 12:38 PM
# @Author  : Runsheng
# @File    : vcf2del.py
"""
parser the vcf file, get a bed like list for deletion
"""
import argparse
import os
import sys
from primervcf.vcfparser import get_deletion_region, bedl2file

parser=argparse.ArgumentParser()
parser.add_argument("-d", "--wkdir", default=None,
                    help="The dir path contain the file, if not given, use the current dir")
# input and out
parser.add_argument("-f", "--file",
                    help="the vcf file used as input")
parser.add_argument("-o", "--out",default="del.out",
                    help="bed file output")
# len cutoff
parser.add_argument("-l", "--length",default=10,
                    help="the min length for a deletion to be used, default is 10")

args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
wkdir=os.getcwd() if args.wkdir is None else args.wkdir

bed_l=get_deletion_region(args.file, key="del", len_cutoff=args.length)
bedl2file(bed_l, filename=args.out)

