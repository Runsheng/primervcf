#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 22/12/2022 9:40 AM
# @Author  : Runsheng
# @File    : vcfparser.py
"""
the vcf related parsers, to convert the vcf to 100% deletion bed
"""
import vcf


def clean_cigar(cigar_string):
    """
    "15b'M'4b'D'7b'M'" to "15M4D7M"
    :param cigar_string:
    :return:
    """
    # remove b
    cigar_1=cigar_string.replace("b", "")
    # remove '
    cigar_2=cigar_1.replace("'", "")
    return cigar_2


def parse_cigar(cigar_string):
    """
    for python3, the return value would be like "15b'M'4b'D'7b'M'"
    need to use clean_cigar to clean the b and '
    :param cigar_string:
    :return:
    """
    #######
    BAM_MATCH = 0
    BAM_INS = 1
    BAM_DEL = 2
    BAM_SOFTCLIP = 4

    CIGAR_OPS = {'M' : BAM_MATCH, 'I' : BAM_INS, 'D' : BAM_DEL, 'S' : BAM_SOFTCLIP}
    ##########
    i = 0
    prev_i = 0
    cigar = []

    cigar_string=clean_cigar(cigar_string)

    while i < len(cigar_string):
        if cigar_string[i] in CIGAR_OPS:
            cigar.append((CIGAR_OPS[cigar_string[i]], int(cigar_string[prev_i:i])))
            prev_i = i + 1
        i += 1
    return cigar


def cigar_getlen(cigar):
    """
    input cigar is a cigar number tuple
    return: match
    """
    M = 0
    for c in cigar:
        if c[0] == 0:
            M += c[1]
    return M


def cigar_getfree3(cigar):
    free3 = 0
    if cigar[-1][0] == 4:
        free3 = cigar[-1][1]

    return free3

### IO functions
def get_deletion_region(vcf_file, key="del", len_cutoff=10):
    """
    get the list for indels: deletion and insertion
    :param vcf_file:
    :param key:
    :return: bed_l: [chro, start, end]
    """

    bed_l=[]
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    for record in vcf_reader:
        #print(record)
        #print(record.var_subtype)
        #print(record.aaf)
        #print(int(record.aaf[0])==1)
        if record.var_subtype.lower()==key and int(record.aaf[0])==1:  # choose only the homo deletion
            alt_seq=record.ALT[0]
            del_len=len(record.REF)-len(alt_seq)
            if del_len>=len_cutoff and "n" not in record.REF: # no gap in the reference
                del_start=record.affected_start+len(alt_seq)-1
                del_end=record.affected_end
                bed_l.append((record.CHROM, del_start,del_end))
    return bed_l


def bedl2file(bed_l, filename):
    with open(filename, "w") as fw:
        for i in bed_l:
            i_str=[str(x) for x in i]
            fw.write("\t".join(i_str))
            fw.write("\n")


def file2bedl(filename):
    """
    chro, start, end, 3col
    :param filename:
    :return:
    """
    bed_l=[]
    with open(filename, "r") as f:
        for line in f.readlines():
            line_str=line.strip().split("\t")
            chro, start, end=line_str
            bed_one=[chro, int(start), int(end)]
            bed_l.append(bed_one)
    return bed_l
