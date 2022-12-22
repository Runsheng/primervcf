#!/usr/bin/env python
#-*- coding: utf-8 -*-
# @Time    : 12/12/2022 12:23 PM
# @Author  : Runsheng
# @File    : primerdesign_vcf.py.py

"""
The main script used to run cmd orders
"""

import argparse
import os
import sys

#currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
#parentdir = os.path.dirname(os.path.dirname(currentdir))
#sys.path.insert(0,parentdir)

from primervcf.primer_indel import  primer3_general_settings
from primerdiffer.utils import fasta2dic,dic2dic

# self import
from primervcf.primer_indel import flow_walk_deletion
from primervcf.vcfparser import file2bedl

parser=argparse.ArgumentParser()
parser.add_argument("-d", "--wkdir", default=None,
                    help="The dir path contain the file, if not given, use the current dir")
# input
parser.add_argument("-g", "--genome",
                    help="the fasta file used to design primer and check specificity")
parser.add_argument("-b", "--bedfile",default=None,
                    help="the bed file containing the deletion region interval")

# primer cutoff
parser.add_argument("--alignlen", type=int, default=16,
                    help="the cutoff of primer min align length as a right hit, default is 16")
parser.add_argument("--free3len",  type=int, default=2,
                    help="the cutoff of primer 3' align length as a right hit, default is 2")

# run parameter, how dense the primer should be
parser.add_argument("-n","--primernumber", type=int, default=5,
                    help="the primer designed for each region, default is 5, do not have much impact for primer design")

parser.add_argument("--debug",  type=str, default="no",
                    help='open debug mode or not, default is no')

# output parameter
parser.add_argument("--prefix", default="primers",
                    help="prefix of output file, default is primers")

# add user settings
parser.add_argument("--primer3config", default=None,
                    help="the config file for the primer3 ")

# default handler
args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
wkdir=os.getcwd() if args.wkdir is None else args.wkdir
debug=False if args.debug=="no" else True

os.chdir(wkdir)

# deal with config
primer3_setting = primer3_general_settings
if args.primer3config is None:
    pass
else: # update the primer3 setting from the given parameters
    with open(args.primer3config, "r") as f:
        primer3_user_setting = eval(f.read()) # read the config as a dict
        primer3_setting.update(primer3_user_setting)

ref_dict=dic2dic(fasta2dic(args.genome))
bed_l=file2bedl(args.bedfile)
#print(bed_l)

flow_walk_deletion(
                   ref_dict=ref_dict,
                   db=args.genome,
                   deletion_bedlist=bed_l,
                   prefix=args.prefix,
                   primer_number=args.primernumber,
                   debugmod=debug,
                   cutoff_alignlength=args.alignlen,
                   cutoff_free3=args.free3len,
                   flank=1000,
                   primer3_general_settings=primer3_setting)
