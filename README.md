# NeoSV
A computational workflow to identify neoantigens from structural variations (SVs)

# Installation
### Prerequisites
* Please make sure your Python version is greater than 3.7.
* NeoSV need netMHCpan to predict the binding affinity between neoantigens and MHC molecules. A detailed guidance for netMHCpan is available here.<br>

### Dowload
  1. via PyPI (recommended): `pip install NeoSV`<br>
  2. via Anaconda: `conda install NeoSV`<br>
  3. or you can download the source code and compile it manually by typing: `python setup.py install`. If you do this way, you also need to install some dependencies, including: biopython, pyensembl.


# User guide
### Input
NeoSV requires 3 types of inputs:
* **_Mutation:_** a file in VCF format which lists all SVs you want to analyze. Please see sv.vcf as a template.
* **_HLA:_** a file listing the HLA types line-by-line. Usually this includes six HLA alleles for a patient. HLA should be in 4 digit format like: HLA-A*02:01. Please see hla.txt as a template.
* **_Reference:_** NeoSV utilizes pyensembl for SV annotation, thus a reference for pyensembl is needed. There are 3 ways to prepare it: <br>
  - **Pre-download by pyensembl (recommended):** When you install NeoSV using pip or conda, pyensembl will be automatically installed as well. Then you can download the reference:<br>
    ```
    export PYENSEMBL_CACHE_DIR=/custom/cache/dir # specify the location for storing reference
    pyensembl install --release <list of Ensembl release numbers> --species <species-name> # download, for hg19 please use release 75, for hg38 please used release 96
    ```
  - **Automatically download by NeoSV:** If NeoSV did not detect a valid reference in , it will automatically download one to that folder. However, please make sure that your server/computer can connect to the internet, because most high performance computing nodes are disconnected.
  - **Prepare the reference file manually:** This would be useful if your data is not from human or mouse. Then you need to prepare the reference by yourself. A FASTA file and a GTF file will be enough. For more details please see the guidance in pyensembl.
### Run
