# NeoSV
<img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/NeoSV">  <a href="https://pypi.python.org/pypi/NeoSV/"> <img src="https://img.shields.io/pypi/v/NeoSV.svg?maxAge=1000" alt="PyPI" /> </a>  <img alt="GitHub" src="https://img.shields.io/github/license/ysbioinfo/NeoSV">  <img alt="GitHub all releases" src="https://img.shields.io/github/downloads/ysbioinfo/NeoSV/total">

A computational workflow to identify **Neo**antigens from **S**tructural **V**ariations.

# New in NeoSV v0.0.3
This release involves 3 updates:
* MHCflurry is enabled as an alternative algorithm for predicting MHC binding. Users can set up MHCflurry following [this guidance](https://github.com/openvax/mhcflurry) and feed the execution path of mhcflurry-predict to `--mhcflurry-path`.
* A transcript list (in Ensembl ID) is allowed as an input `--transcript-list`. NeoSV will first use the transcripts in this list when there are multiple isoforms, instead of using the longest one.
* Accept lower cases of [ATCG] in the REF and ALT field of input VCF.

# Background
Neoantigens are considered as ideal targets for immunotherapies because they are tumor-specifc and not subject to immune tolerance. Previous studies have been focused on single nucleotide variation (SNV) and insertion-and-deletion (indel), with the neoantigens from structural variation (SV) less characterized.

We developed a Python package-NeoSV-to **_annotate_** the effect of SVs on protein and **_predict_** potential neoantigens created by SVs. We have successfully applied NeoSV to thousands of tumor genomes from Pan Cancer Analysis of Whole Genomes (PCAWG) and constructed a comprehensive repertoire of SV-derived neoantigens. For more detailed information, please see our paper. 

# Install
### Prerequisites
* [Python (>3.6)](https://www.python.org/downloads/). NeoSV should work well with all versions of Python3, but has been only tested on Python > 3.6
* [NetMHCpan (>4.0)](https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1). After you sign up and get the link for downloading, there will be a accompanied guidance on how to configure netMHCpan.
* [MHCflurry (2.0)](https://github.com/openvax/mhcflurry) (optional). Anaconda may help on installing because MHCflurry requires tensorflow.
### Download
* PyPI: if you already have python and pip, you can directly install NeoSV via `pip install neosv`<br>
* Source code: we noted that sometimes pip will not install the binary file neosv, is such case you can download the package and install it using `python setup.py install`. Please remember to install [biopython](https://biopython.org/) and [pyensembl](https://github.com/openvax/pyensembl) using pip before installation.


# Usage
### Input
NeoSV requires 3 types of inputs:
* **Mutation file:** a file in [VCF format](https://samtools.github.io/hts-specs/VCFv4.2.pdf) which lists all SVs you want to analyze. Please see [test.sv.vcf](https://github.com/ysbioinfo/NeoSV/blob/main/test.sv.vcf) as a template.
* **HLA file:** a file listing the HLA types line-by-line. Usually this includes six HLA alleles for a patient. HLA should be in 4 digit format like: HLA-A*02:01. Please see [test.hla.txt](https://github.com/ysbioinfo/NeoSV/blob/main/test.hla.txt) as a template.
* **Reference file:** NeoSV utilizes pyensembl for SV annotation, thus a reference for pyensembl is needed. There are 3 ways to prepare it: <br>

  - _Pre-download by pyensembl (recommended):_ When you install NeoSV using pip or conda, pyensembl will be automatically installed as well. Then you can download the reference:<br>
    ```
    export PYENSEMBL_CACHE_DIR=/custom/cache/dir # specify the location for storing reference
    pyensembl install --release <list of Ensembl release numbers> --species <species-name> # download, for hg19 please use release 75, for hg38 please used release 96
    ```
  - _Automatically download by NeoSV:_ If NeoSV did not detect a valid reference in `--pyensembl-cache-dir`, it will automatically download one to that folder. However, please make sure that your server/computer can connect to the internet, because most high performance computing nodes are disconnected.
  - _Prepare the reference file manually:_ This would be useful if your data is not from human or mouse. Then you need to prepare the reference by yourself. A FASTA file and a GTF file will be enough. For more details please see the [guidance](https://github.com/openvax/pyensembl#non-ensembl-data).
### Run
* Quick start: suppose you have a mutation file named `test.sv.vcf`, a HLA file named `test.hla.txt`. Your pyensembl reference is human sapiens release 75 and located at `/pyensembl/`, then a typical NeoSV command is:
  ```
  neosv -vf test.sv.vcf -hf test.hla.txt -np /path/to/netmhcpan -o test -p test -r 75    
  ```
* Below is detailed description for each parameter:
  | &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  &nbsp;&nbsp;&nbsp;&nbsp;  &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; Argument &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  &nbsp;&nbsp;&nbsp;&nbsp;  &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;| Description |
  | :------------ | --- |
  | `-h`, `--help` | show the help message |
  | `-vf`, `--vcf-file` | Structural variants in VCF format. |
  | `-hf`, `--hla-file` | HLA alleles (resolution: 4 digit), with one allele per line. |
  | `-np`, `--netmhc-path` | Absolute path to the NetMHCpan execution file, please skip this argument if NetMHCpan has been added to your PATH. |
  | `-mp`, `--mhcflurry-path` | Absolute path to the MHCflurry execution file. This argument is optional if `-np` has been specified. |
  | `-o`, `--out` | Folder for all result files. A new folder will be created if it does not exist. |
  | `-p`, `--prefix` | This prefix will be added to all output files. |
  | `-r`, `--release` | The release of Ensembl to use. Valid release versions can be found here. Ensembl release corresponding to hg19/GRCh37, hg38/GRCh38 are 75, 95. If your data is from other species, you need to download a GTF file and a cDNA file from [Ensembl](ftp://ftp.ensembl.org/pub) and specify them using -gf and -cf |
  | `-gf`, `--gtf-file` | GTF file for the reference, only needed  |
  | `-cf`, `--cdna-file` | cDNA file for the reference |
  | `-pd`, `--pyemsembl-cache-dir` | Directory for Pyensembl cache files. If not specified, the platform-specific cache folder will be used |
  | `-l`, `--epitope-lengths` | Lengths of neoepitopes to predict MHC binding. Default: 8-11. |
  | `-ic`, `--ic50-cutoff` | Filter neoepitopes with IC50 (nM) above this value. Default: 500. |
  | `-rc`, `--rank-cutoff` | Filter neoepitopes with rank above this value. Default: 2. |
  | `-ct`, `--complete-transcript` | Only complete transcripts will be considered for SV annotation. Default: True. |
  | `-t`, `--transcript-list` | A list of transcripts (in Ensembl ID) with the highest priority. Default: None. |
  | `--anno-only` | Only annotate SV without predicting neoantigens.If this argument is added, --hla-file is not required, and you will only get the annotation result. |

### Output
Several files will be generated in the output directory, you may have interest in the files suffixed by ***neoantigen.filtered.txt*** and ***anno.filtered.txt***

* ***{prefix}.neoantigen.filtered.txt*** stores all information of the candidate neantigens:

  | Column index | Column name | Content |
  | :------------: | --- | --- |
  | 1 | chrom1 | Chromosome of the 1st breakpoint |
  | 2 | pos1 | Genommic position of the 1st breakpoint |
  | 3 | gene1 | Gene name of the 1st breakpoint |
  | 4 | transcript_id1 | Ensembl transcript ID of the 1st breakpoint |
  | 5 | chrom2 | Chromosome of the 2nd breakpoint |
  | 6 | pos2 | Genommic position of the 2nd breakpoint |
  | 7 | gene2 | Gene name of the 2nd breakpoint |
  | 8 | transcript_id2 | Ensembl transcript ID of the 2nd breakpoint |
  | 9 | svpattern | |
  | 10 | svtype | SV types according to the orientation of junction read. Values: _DUP_, _DEL_, _TRA_, _t2tINV_, or _h2hINV_. |
  | 11 | frameshift | The effect on open reading frame. Values: _In-frame_, _Stop-gain_, _Stop-loss_, _Start-loss_. |
  | 12 | neoantigen | Amino acid sequence of the neoantigen |
  | 13 | allele | HLA allele that binds to the neoantigen |
  | 14 | affinity | Binding affinity (nM) provided by netMHCpan |
  | 15 | rank | Rank of the binding provided by netMHCpan |

* ***{prefix}.anno.filtered.txt*** stores all annotations of the SVs:

  | Column index | Column name | Content |
  | :------------: | --- | --- |
  | 1 | chrom1 | Chromosome of the 1st breakpoint. |
  | 2 | pos1 | Genommic position of the 1st breakpoint. |
  | 3 | function1 | Location of the 1st breadpoint relative to a gene. Values: _Intergenic_, _Intron_, _Exon_. |
  | 4 | gene1 | Gene name of the 1st breakpoint. |
  | 5 | transcript_id1 | Ensembl transcript ID of the 1st breakpoint |
  | 6 | strand1 | Coding strand of the 1st gene. Values: _+_, _-_, _None (if intergenic)_ |
  | 7 | transcript_retain1 | The part being retained of transcript, I/i indicates intron, E/e indicates exon. Upper case means an intact exon/intron, while lower case means the exon/intron is truncated by this SV |
  | 8 | chrom2 | Chromosome of the 2nd breakpoint |
  | 9 | pos2 | Genommic position of the 2nd breakpoint |
  | 10 | function2 | Location of the 1st breadpoint relative to a gene. Values: _Intergenic_, _Intron_, _Exon_. |
  | 11 | gene2 | Gene name of the 2nd breakpoint |
  | 12 | transcript_id2 | Ensembl transcript ID of the 2nd breakpoint |
  | 13 | strand2 | Coding strand of the 2nd gene. Values: _+_, _-_, _None (if intergenic)_ |
  | 14 | transcript_retain2 | |
  | 15 | svpattern | |
  | 16 | svtype | SV types according to the orientation of junction read. Values: _DUP_, _DEL_, _TRA_, _t2tINV_, or _h2hINV_. |
  | 17 | fusion | Whether this SV can lead to a functional gene fusion. It should be noted that the fusion is not restricted to two-gene fusion. |

* ***{prefix}.net.in.txt*** stores the peptides fed to netMHCpan.

* ***{prefix}.net.out.txt*** stores the raw output from netMHCpan.

# License
NeoSV is licensed under the terms of MIT license.
