#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 5/12/2022 8:29 AM
# @Author  : Runsheng
# @File    : map2vcf.py


"""
The pipeline wrapper to get a vcf file from fastq mapping, need:
1. bwa mem
2. samtools
3. vcftools
"""
import os

from primervcf.utils import myexe

def wrapper_bwamem(index, read_list, prefix="default", core=32 , wkdir=None):
    """
    general wrapper for bwa mem to be used in python cmd line
    para: index, the bwa index of reference fasta file
    para: read_list, the fastq file names in a python list
    para: prefix, the prefix for the output bam file
    para: core, the cores used for mapping and sorting

    return: None
    """
    prefix = read_list[0].split("_")[0] if prefix == "default" else prefix
    read_str = " ".join(read_list)

    if wkdir is None:
        wkdir=os.getcwd()
    os.chdir(wkdir)

    myexe("bwa mem -t {core} {index} {read_str} > {prefix}.sam".format(core=core, index=index, read_str=read_str,
                                                                       prefix=prefix))
    myexe("samtools view -bS {prefix}.sam >{prefix}.bam".format(prefix=prefix))
    myexe("samtools sort -@ {core} {prefix}.bam {prefix}_s".format(core=core, prefix=prefix))
    myexe("samtools index {prefix}_s.bam".format(prefix=prefix))
    myexe("rm {prefix}.sam {prefix}.bam".format(prefix=prefix))

    print ("++++++++done++++++++++")
    return prefix+"_s.bam"


def wrapper_bam2vcf(ref, bamfile, prefix="default", wkdir=None):
    """
    para: ref, the reference fasta file
    para: bamfile: a single sorted bam file
    para: prefix, the prefix for the output vcf file

    return: None
    """
    prefix = bamfile.split(".")[0].split("_")[0] if prefix == "default" else prefix

    if wkdir is None:
        wkdir=os.getcwd()
    os.chdir(wkdir)
    cmd_mpileup = ("bcftools mpileup -Ob -o {prefix}.bcf QR25.bcf -f {ref} {bamfile}"
                   .format(ref=ref, bamfile=bamfile, prefix=prefix))
    cmd_tovcf = ("bcftools call -vmO z -o {prefix}.vcf {prefix}.bcf".
                 format(prefix=prefix))

    print(cmd_mpileup)
    myexe(cmd_mpileup)
    print(cmd_tovcf)
    myexe(cmd_tovcf)


def flow_fastq2vcf(ref,read, prefix="raw", core=16 ,wkdir=None):
    """
    read=["/home/zhaolab1/data/dnaseq/trim/AF16_F_P.fq.gz", "/home/zhaolab1/data/dnaseq/trim/AF16_R_P.fq.gz"]
    prefix="AF16"
    ref="/home/zhaolab1/data/dnaseq/refindex/cb4.fa"
    wrapper_bwamem(ref, read, prefix)
    bamfile="AF16_s.bam"
    ref="/home/zhaolab1/data/dnaseq/refindex/cb4.fa"
    wrapper_bam2vcf(ref, bamfile)

    # use SRR1793006 as example, QR25 fastq

    :return:
    """
    if wkdir is None:
        wkdir=os.getcwd()
    os.chdir(wkdir)
    prefix = read[0].split("_")[0] if prefix == "raw" else prefix
    bamfile=wrapper_bwamem(ref, read, prefix, core=core)
    wrapper_bam2vcf(ref, bamfile)
