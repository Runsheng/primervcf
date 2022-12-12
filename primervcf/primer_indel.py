# !/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 28/10/2022 12:47 PM
# @Author  : Runsheng
# @File    : primer_indel.py

"""
use the long indel information to design hyplotype specific primers
"""

from primerdiffer.utils import reverse_complement, chr_select
from primerdiffer.primer_check import insilicon_pcr
from primerdiffer.general_settings import primer3_general_settings

#third part import
import primer3 # primer3-py
from pyssw.ssw_wrap import Aligner
import vcf #pyvcf, use pyVCF3 instead of original pyvcf


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


def has_local_hit(primer, seq, cutoff_alignlength, cutoff_free3, debugmod=False):
    aligner = Aligner(seq, report_cigar=True)
    aln = aligner.align(primer)
    cigar = parse_cigar(aln.cigar_string)

    M = cigar_getlen(cigar)
    free3 = cigar_getfree3(cigar)

    if debugmod == True:
        print(aln.score, aln.ref_begin, aln.ref_end, aln.query_begin, aln.query_end, aln.cigar_string)

    if M <= cutoff_alignlength or free3 >= cutoff_free3:
        if debugmod:
            print("M", M)
            print("free3",free3)
        return False
    else:
        return True


def has_primer_local_hit(primer_left, primer_right, seq,
                         cutoff_alignlength=16, cutoff_free3=2,  debugmod=False):
    """
    para: the left and right primers
    return: a bed-like tuple-list

    need: parser_cigar, cigar_getlen, cigar_getfree3
    """
    # F or R primer do not have hit, then it works
    if has_local_hit(primer_left, seq, cutoff_alignlength, cutoff_free3, debugmod) == True and \
            has_local_hit(reverse_complement(primer_right), seq, cutoff_alignlength, cutoff_free3, debugmod) == True:
        return True
    else:
        return False


def my_design_primer_del(name, seq, primer3_general_settings, flank=1000):
    # the flanking is used to ensure the primer would cover the deletion region
    seq_args = {'SEQUENCE_ID': name,
                'SEQUENCE_TEMPLATE': seq,
                'SEQUENCE_OVERLAP_JUNCTION_LIST': flank,
                'PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION': 6
                }
    # IMPORTANT: do not misuse the seq_args and the general_args!
    myprimer = primer3.bindings.designPrimers(seq_args, primer3_general_settings)
    return myprimer


############ IO functions
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
        if record.var_subtype==key and record.aaf[0]==1:  # choose only the homo deletion
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
############


def del_primer_check(del_single, db, ref_dict, primer_number=5, debugmod=False,
                     primer3_general_settings=primer3_general_settings, flank=1000,
                     cutoff_alignlength=16,cutoff_free3=2):
    '''
    print del_primer_check(del_QR25[1], db="/home/zhaolab1/data/dnaseq/refindex/cb4",ref_dict=cb4, debugmod=True)

    del_single, use this to design a primer, have bed like format as chro, start, end
    db is blastdb

    call:
    insilicon_pcr: have only one match in current genome
    check_deletion_local： no local match for the deleted sequence
    return 0 means the primer can be found
    '''
    chro, start, end = del_single
    del_len = end - start
    name, seq = chr_select(ref_dict, chro, start - 1000, end + 1000)  # the origin sequence
    seq_new = seq[:1000] + seq[(1000 + del_len):]  # the seq with deletion inside
    myprimer = my_design_primer_del(name=name, seq=seq,
                                    primer3_general_settings=primer3_general_settings,
                                    flank=flank)

    product_l = []

    for i in range(0, primer_number):
        try:
            left = myprimer['PRIMER_LEFT_' + str(i) + '_SEQUENCE']
            right = myprimer['PRIMER_RIGHT_' + str(i) + '_SEQUENCE']
            product_size = myprimer['PRIMER_PAIR_' + str(i) + '_PRODUCT_SIZE']
            product_l = insilicon_pcr(left, right, db)
        except KeyError:
            pass

        if debugmod:
            print("The %d primer :" % i)
            print("in silicon PCR product number:", len(product_l))
            print ("L:{}, R:{}, rcL:{}, rcV:{}".format( left, right, reverse_complement(left), reverse_complement(right)) )
            print ("has_primer_local_hist in delseq, should be False", has_primer_local_hit(left, right, seq_new, debugmod=debugmod))
            print ("has_primer_local_hit in rawseq, should be True", has_primer_local_hit(left, right, seq, debugmod=debugmod))

        if len(product_l) == 1 and has_primer_local_hit(left, right, seq_new, cutoff_alignlength=cutoff_alignlength, cutoff_free3=cutoff_free3,  debugmod=debugmod)==False:
            # print ("pass")
            # with open((name+".txt"), "w") as fw:
            # fw.write("%s\t%s\t%s\t%d" % (name,left,right,product_size))
            return (name, left, right, product_size)
        else:
            pass
    return 0


def flow_walk_deletion(ref_dict, db, deletion_bedlist, prefix="primers",
                       primer_number=5, debugmod=False,
                       primer3_general_settings=primer3_general_settings, flank=1000,
                       cutoff_alignlength=16, cutoff_free3=2
                       ):
    """
    primer_dict=walk_deletion(ref_dict=cb4,db="/home/zhaolab1/data/dnaseq/refindex/cb4",deletion_bedlist=del_QR25)

    :param ref_dict:
    :param db:
    :param deletion_bedlist:
    :param prefix:
    :return:
    """
    primer_dict = {}
    f_out = open(prefix + ".txt", "w")

    left, right, product_size = (0, 0, 0)
    for del_single in deletion_bedlist:
        primer_status = del_primer_check(del_single, db, ref_dict,
                     primer_number=primer_number, debugmod=debugmod,
                     primer3_general_settings=primer3_general_settings, flank=flank,
                     cutoff_alignlength=cutoff_alignlength,cutoff_free3=cutoff_free3)
        if primer_status == 0:
            print("no primer get for", del_single)
        else:
            name, left, right, product_size = primer_status
            print("Get primer in %s!" % name)
            primer_dict[name] = (left, right, product_size)
            f_out.write(name + "\t" + left + "\t" + right + "\t" + str(product_size) + "\n")
    f_out.close()
    return primer_dict
