#!/usr/bin/env python
#-*- coding: utf-8 -*-
# @Time    : 12/12/2022 1:25 PM
# @Author  : Runsheng
# @File    : fq2vcf.py.py

"""
mapping fastq to reference with bwa mem, and call the vcf
"""
import argparse
import os
import sys

from primervcf.map2vcf import flow_fastq2vcf

parser=argparse.ArgumentParser()
parser.add_argument("-d", "--wkdir", default=None,
                    help="The dir path contain the file, if not given, use the current dir")
# input and out
parser.add_argument("-f", "--file",
                    help="the fastq file used as input, name seperated by , "
                         "example for pair-end like SRR1793006_1.fastq,SRR1793006_2.fastq")
parser.add_argument("-g", "--genome",
                    help="the genome file used for mapping")
parser.add_argument("-p", "--prefix",default="raw",
                    help="the outvcf filename prefix, default is using the file prefix of fastq")

# performace
parser.add_argument("-c", "--core",
                    help="the core used for bwa mem mapping and samtools sort")


args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
wkdir=os.getcwd() if args.wkdir is None else args.wkdir

os.chdir(wkdir)
name_str=args.file
read_list=name_str.strip().split(",")

flow_fastq2vcf(ref=args.genome,
               read=read_list,
               prefix=args.prefix,
               core=args.core,
               wkdir=args.wkdir
               )