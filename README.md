# NeoSV
<img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/NeoSV">  <a href="https://pypi.python.org/pypi/NeoSV/"> <img src="https://img.shields.io/pypi/v/NeoSV.svg?maxAge=1000" alt="PyPI" /> </a>  <img alt="GitHub" src="https://img.shields.io/github/license/ysbioinfo/NeoSV">  

A computational workflow to identify **Neo**antigens from **S**tructural **V**ariations.

# Background
Neoantigens are considered as ideal targets for immunotherapies because they are tumor-specifc and not subject to immune tolerance. Previous studies have been focused on single nucleotide variation (SNV) and insertion-and-deletion (indel), with the neoantigens from structural variation (SV) less characterized.

We developed a Python package-NeoSV-to **_annotate_** the effect of SVs on protein and **_predict_** potential neoantigens created by SVs. We have successfully applied NeoSV to thousands of tumor genomes from Pan Cancer Analysis of Whole Genomes (PCAWG) and constructed a comprehensive repertoire of SV-derived neoantigens. For more detailed information, please see our paper. 

# Install
### Prerequisites
* [Python (>3.6)](https://www.python.org/downloads/). NeoSV should work well with all versions of Python3, but has been only tested on Python > 3.6
* [netMHCpan (>4.0)](https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1). After you sign up and get the link for downloading, there will be a accompanied guidance on how to configure netMHCpan.
### Download
* PyPI (recommended): if you already have python and pip, you can directly install NeoSV via `pip install NeoSV`<br>
* From source: download the code from github and type `python setup.py install`. If you do this way, you need install [biopython](https://biopython.org/), [pyensembl](https://github.com/openvax/pyensembl) by yourself.


# User guide
### Input
NeoSV requires 3 types of inputs:
* **Mutation file:** a file in VCF format which lists all SVs you want to analyze. Please see sv.vcf as a template.
* **HLA file:** a file listing the HLA types line-by-line. Usually this includes six HLA alleles for a patient. HLA should be in 4 digit format like: HLA-A*02:01. Please see hla.txt as a template.
* **Reference file:** NeoSV utilizes pyensembl for SV annotation, thus a reference for pyensembl is needed. There are 3 ways to prepare it: <br>
  - **Pre-download by pyensembl (recommended):** When you install NeoSV using pip or conda, pyensembl will be automatically installed as well. Then you can download the reference:<br>
    ```
    export PYENSEMBL_CACHE_DIR=/custom/cache/dir # specify the location for storing reference
    pyensembl install --release <list of Ensembl release numbers> --species <species-name> # download, for hg19 please use release 75, for hg38 please used release 96
    ```
  - **Automatically download by NeoSV:** If NeoSV did not detect a valid reference in , it will automatically download one to that folder. However, please make sure that your server/computer can connect to the internet, because most high performance computing nodes are disconnected.
  - **Prepare the reference file manually:** This would be useful if your data is not from human or mouse. Then you need to prepare the reference by yourself. A FASTA file and a GTF file will be enough. For more details please see the guidance in pyensembl.
### Run
Suppose you have a mutation file named `test.sv.vcf`, a HLA file named `test.hla.txt`. Your pyensembl reference is human sapiens release 75 and located at `/pyensembl/`, then a typical NeoSV command is:
```
python neosv.py -vf test.sv.vcf -hf test.hla.txt -np /path/to/netmhcpan -o test -p test -r 75    
```
Below is detailed description for each parameter:
| Parameter | Description |
| :------------: | --- |
| -h | show the help message |
| -vf | File for structural variants in VCF format |
| -hf | File for HLA alleles with one allele per line |
| -np | Absolute path of the netMHCpan execution file |
| -o | Output directory. A new folder will be created if it doesn't exist |
| -p | The prefix added to all output files |
| -r | Which reference (ENSEMBL release) you want to use. Valid values are54-94(homo sapiens release number) or custom. Ensembl release corresponding to hg19/GRCh37, hg38/GRCh38 are 54, 75, 95.If your data are from other species(custom), please download the gtf file and the cdna file from ENSEMBL we ftp://ftp.ensembl.org/pub and specify them using --gtf-file and --cdna-file |
| -gf | Path of GTF file for the assembly |
| -cf | Path of CDNA file for the assembly |
| -pd | Directory for Pyensembl cache files. If not specified, the paltform-specificCache folder will be used |
| -l | Lengths of neoepitopes to predict MHC binding |
| -ic | Filter neoepitopes with IC50 (nM) above this value |
| -rc | Filter neoepitopes with rank above this value |
| -ct | Only complete transcripts will be considered for SV annotation |
| --anno-only | Whether to only annotate SV without predicting neoantigens |

### Output

# Citation
