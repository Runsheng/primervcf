#  **primervcf** package
![PyPI](https://img.shields.io/pypi/v/primervcf?color=green)

## Installation:
The package worked with python version >=3.4.
Only tested in linux x64 system.

python package:
- primer3-py>=0.6.1
- pyssw >=0.1.4
- primerdiffer >=0.1.4
- biopython>=1.7.8
- pyssw>=0.1.4

other program in your $PATH:
- ncbi-blast

example code to install the packages under python3 with pip and conda
```bash
pip install primervcf # will also install all other python dependencies 
conda install -c bioconda blast # install ncbi blast, which is not included in pip installation
```

## Walkthrough
Design primers based on the VCF file provided, parser the deletion, and force primers to overalp the deletion 
- Check all long indel regions (>=indel_cutoff, default is 10) to design primers.
- Check local region to make sure the vcf is not 
- Use the whole genome (--genome1) as db to make a specificity check, to ensure the primer can only amplify one region.
- The dis-similarity between genome and the strain could be individual/strain/population level (<=1%).   

Use _C.briggsae_ QR25 strain as example. The fasta files and the vcf can be downloaded separately 
from [cb5.fa](https://github.com/Runsheng/cbgenome/releases/download/cb5pre_cn3pre/cb5.fa.gz) and 
from [QR25.vcf]()


Quick example
```bash
vcf2del.py -f QR25.vcf -o QR25_del.bed
primerdesign_vcf.py -g cb5.fa -b QR25_del.bed --prefix QR25
```


Explanations:

```bash
usage: vcf2del.py [-h] [-d WKDIR] [-f FILE] [-o OUT] [-l LENGTH]
optional arguments:
  -h, --help            show this help message and exit
  -d WKDIR, --wkdir WKDIR
                        The dir path contain the file, if not given, use the current dir
  -f FILE, --file FILE  the vcf file used as input
  -o OUT, --out OUT     bed file output
  -l LENGTH, --length LENGTH
                        the min length for a deletion to be used, default is 10
# run example
vcf2del.py -f QR25.vcf -o QR25_del.bed

# check the result in file "QR25_del.bed"
head QR25_del.bed
ChrI    144730  144744
ChrI    423135  423167
ChrI    461207  461222
ChrI    465639  465661
ChrI    492790  492800

# run primerdesign_vcf.py to output primer
usage: primerdesign_vcf.py [-h] [-d WKDIR] [-g GENOME] [-b BEDFILE] [--alignlen ALIGNLEN] [--free3len FREE3LEN]
                           [-n PRIMERNUMBER] [--debug DEBUG] [--prefix PREFIX] [--primer3config PRIMER3CONFIG]

optional arguments:
  -h, --help            show this help message and exit
  -d WKDIR, --wkdir WKDIR
                        The dir path contain the file, if not given, use the current dir
  -g GENOME, --genome GENOME
                        the fasta file used to design primer and check specificity
  -b BEDFILE, --bedfile BEDFILE
                        the bed file containing the deletion region interval
  --alignlen ALIGNLEN   the cutoff of primer min align length as a right hit, default is 16
  --free3len FREE3LEN   the cutoff of primer 3' align length as a right hit, default is 2
  -n PRIMERNUMBER, --primernumber PRIMERNUMBER
                        the primer designed for each region, default is 5, do not have much impact for primer design
  --debug DEBUG         open debug mode or not, default is no
  --prefix PREFIX       prefix of output file, default is primers
  --primer3config PRIMER3CONFIG
                        the config file for the primer3

# run example, chose only the first 10 deletion to design primers
head QR_25_del.bed > QR25_delhead.bed
primerdesign_vcf.py -g cb5.fa -b QR25_delhead.bed --prefix QR25
```



## Control parameters: design primers with given parameters
The default primer design parameter is described in [general_setting.py](https://github.com/Runsheng/primerdiffer/blob/master/primerdiffer/general_settings.py):
```python
primer3_general_settings =  {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_PICK_LEFT_PRIMER':1,
        'PRIMER_PICK_RIGHT_PRIMER':1,
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 23,
        'PRIMER_OPT_TM': 57.0,
        'PRIMER_MIN_TM': 46.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_PRODUCT_SIZE_RANGE': [[250, 650]],
        'PRIMER_NUM_RETURN':10,
        'PRIMER_MIN_THREE_PRIME_DISTANCE':10,
        'PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS':0
}
```
The new parameter can be supplemnted as a file with the terms which need to be changed, for example, create a file with name of config.txt
```python
# content of config.txt
# only changed the product size, the others are the same the default
{'PRIMER_PRODUCT_SIZE_RANGE': [[100, 200]]}
```


Other scripts which might be useful if you want to start from fastq mapping
```bash
usage: fq2vcf.py [-h] [-d WKDIR] [-f FILE] [-g GENOME] [-p PREFIX] [-c CORE]

optional arguments:
  -h, --help            show this help message and exit
  -d WKDIR, --wkdir WKDIR
                        The dir path contain the file, if not given, use the current dir
  -f FILE, --file FILE  the fastq file used as input
  -g GENOME, --genome GENOME
                        the genome file used for mapping
  -p PREFIX, --prefix PREFIX
                        the outvcf filename prefix, default is using the file prefix of fastq
  -c CORE, --core CORE  the core used for bwa mem mapping and samtools sort, default is 16

# this script require the bwa mem in your $PATH
# which can be installed by "conda install -c bioconda bwa"
# example
fq2vcf.py -f SRR1793006_1.fastq,SRR1793006_2.fastq -g cb5.fa -c 32
```



    


            