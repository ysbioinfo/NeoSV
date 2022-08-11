# NeoSV
A computational workflow to identify neoantigens from structural variations (SVs)

# Installation
Prerequisite:
* Please make sure your Python version is greater than 3.7.
* NeoSV need netMHCpan to predict the binding affinity between neoantigens and MHC molecules. A detailed guidance for netMHCpan is available here. 
Then you can install NeoSV in three ways:<br>
1. via PyPI (recommended): `pip install NeoSV`<br>

2. via Anaconda: `conda install NeoSV`<br>

3. or you can download the source code and compile it manually by typing: `python setup.py install`. If you do this way, you also need to install some dependencies, including: biopython, pyensembl.


# User guide
## Input
NeoSV requires 3 types of inputs:
* Mutation: a file in VCF format which lists all SVs you want to analyze. Please see sv.vcf as a template.
* HLA: a file listing the HLA types line-by-line. Usually this includes six HLA alleles for a patient. HLA should be in 4 digit format like: HLA-A*02:01. Please see hla.txt as a template.
* Reference: NeoSV utilizes pyensembl for SV annotation, so you need to prepare reference files:

## Run
