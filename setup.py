from setuptools import setup
from primervcf import __version__

# read the description from the readme file
with open('readme.md', "r") as f:
    LONG_DESC = f.read()

setup(
    name='primervcf',
    version=__version__,
    packages=['', "primervcf", "bin"],
    url='https://github.com/Runsheng/primervcf',
    license='GPL-2',
    author='runsheng',
    author_email='runsheng.lee@gmail.com',
    description='primer design for haplotype genotyping using indel information',
    long_description=LONG_DESC,
    long_description_content_type='text/markdown',
    install_requires = ["primer3-py>=0.6.1",
                        "biopython>=1.78",
                        "primerdiffer>=0.1.4",
                        "pyssw>=0.1.4",
                        "PyVCF3>=1.0.3"],
    scripts = ['bin/primerdesign_vcf.py',
               'bin/vcf2del.py',
               'bin/fq2vcf.py'
                ]
)


