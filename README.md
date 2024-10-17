# NeoSV
<img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/NeoSV">  <a href="https://pypi.python.org/pypi/NeoSV/"> <img src="https://img.shields.io/pypi/v/NeoSV.svg?maxAge=1000" alt="PyPI" /> </a>  <img alt="GitHub" src="https://img.shields.io/github/license/ysbioinfo/NeoSV">

A computational workflow to identify **Neo**antigens from **S**tructural **V**ariations.

# New in NeoSV v0.0.4
* Support BEDPE format as input
* Fix bugs related to NetMHCpan 4.1 (NetMHCpan 4.0 will no longer be supported by NeoSV.)
* Add an additional parameter erc, which enable users filter neoantigens by EL (eluted ligand) rank

# Background
Neoantigens are considered as ideal targets for immunotherapies because they are tumor-specifc and not subject to immune tolerance. Previous studies have been focused on single nucleotide variation (SNV) and insertion-and-deletion (indel), with the neoantigens from structural variation (SV) poorly characterized.

We developed a Python package-NeoSV-to **_annotate_** the effect of SVs on protein and **_predict_** potential neoantigens created by SVs. We have successfully applied NeoSV to thousands of tumor genomes from Pan Cancer Analysis of Whole Genomes (PCAWG) and constructed a comprehensive repertoire of SV-derived neoantigens. For more details, please read our paper: 

> Shi, Y., Jing, B. & Xi, R. Comprehensive analysis of neoantigens derived from structural variation across whole genomes from 2528 tumors. Genome Biol 24, 169 (2023)


# Install
### Prerequisites
* [Python (>3.6)](https://www.python.org/downloads/). NeoSV should work well with all versions of Python3, but has been only tested on Python > 3.6
* [NetMHCpan (4.1)](https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1). After you sign up and get the link for downloading, there will be a accompanied guidance on how to configure netMHCpan.
### Download
* PyPI: if you already have python and pip, you can directly install NeoSV via `pip install neosv`<br>
* Source code: we noted that sometimes pip will not install the binary file neosv, is such case you can download the package and install it using `python setup.py install`. Please remember to install [biopython](https://biopython.org/) and [pyensembl](https://github.com/openvax/pyensembl) using pip before installation.


# Usage
### Input
NeoSV requires 3 types of inputs:
* **Variant file:** a file in [VCF format](https://samtools.github.io/hts-specs/VCFv4.2.pdf) or [BEDPE format](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format) which lists all SVs you want to analyze. Template files: [test.sv.vcf](https://github.com/ysbioinfo/NeoSV/blob/main/test.sv.vcf) and [test.sv.bedpe](https://github.com/ysbioinfo/NeoSV/blob/main/test.sv.bedpe)
* **HLA file:** a file listing the HLA alleles line by line. This usually includes six HLA alleles for an individual. HLA should be in 4 digit format like HLA-A*02:01. Template file: [test.hla.txt](https://github.com/ysbioinfo/NeoSV/blob/main/test.hla.txt)
* **Reference file:** NeoSV utilizes pyensembl for SV annotation, thus a reference for pyensembl is needed. There are 3 ways to prepare it: <br>

  - _Pre-download by pyensembl (recommended):_ When you install NeoSV using pip or conda, pyensembl will be automatically installed as well. Then you can download the reference:<br>
    ```
    export PYENSEMBL_CACHE_DIR=/custom/cache/dir # specify the location for storing reference
    pyensembl install --release <list of Ensembl release numbers> --species <species-name> # download, for hg19 please use release 75, for hg38 please used release 96
    ```
  - _Automatically download by NeoSV:_ If NeoSV did not detect a valid reference in `--pyensembl-cache-dir`, it will automatically download one to that folder. Please make sure the internet connection of your system, since some high performance computing nodes have no network.
  - _Prepare the reference file manually:_ This would be useful if your data is not from human or mouse. Then you need to prepare the reference by yourself. A FASTA file and a GTF file will be enough. For more details please see the [guidance](https://github.com/openvax/pyensembl#non-ensembl-data). In addition, you need to confirm the MHC alleles in that species are supported by NetMHCpan. 
### Run
* Quick start: suppose you have a variant file named `test.sv.vcf`, a HLA file named `test.hla.txt`. Your pyensembl reference is human sapiens release 75 and located at `/pyensembl/`, then a typical NeoSV command is:
  ```
  neosv -sf test.sv.vcf -hf test.hla.txt -np /path/to/netmhcpan -o test -p test -r 75    
  ```
* Below is detailed description for each parameter:
  | &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  &nbsp;&nbsp;&nbsp;&nbsp;  &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; Argument &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  &nbsp;&nbsp;&nbsp;&nbsp;  &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;| Description |
  | :------------ | --- |
  | `-h`, `--help` | show the help message |
  | `-sf`, `--sv-file` | Structural variants in VCF or BEDPE format. NeoSV will automatically identify the format according to the file suffix.|
  | `-hf`, `--hla-file` | HLA alleles (resolution: 4 digit), with one allele per line. |
  | `-np`, `--netmhc-path` | Absolute path to the NetMHCpan execution file, please skip this argument if NetMHCpan has been added to your PATH. |
  | `-o`, `--out` | Folder for all result files. A new folder will be created if it does not exist. |
  | `-p`, `--prefix` | This prefix will be added to all output files. |
  | `-r`, `--release` | The release of Ensembl to use. Valid release versions can be found here. Ensembl release for hg19/GRCh37, hg38/GRCh38 are 75, 96. |
  | `-gf`, `--gtf-file` | GTF file for the reference, only needed when you want to prepare the ensembl reference by yourself. |
  | `-cf`, `--cdna-file` | cDNA file for the reference, only needed when you want to prepare the ensembl reference by yourself. |
  | `-pd`, `--pyemsembl-cache-dir` | Directory for Pyensembl cache files. If not specified, the platform-specific cache folder will be used |
  | `-l`, `--epitope-lengths` | Lengths of neoepitopes to predict MHC binding. Default: 8-11. |
  | `-ic`, `--ic50-cutoff` | Filter neoepitopes with IC50 (nM) above this value. Default: 500. |
  | `-brc`, `--ba-rank-cutoff` | Filter neoepitopes with BA-rank above this value. Default: 2. |
  | `-erc`, `--el-rank-cutoff` | Filter neoepitopes with EL-rank above this value. Default: 2. |
  | `-ct`, `--complete-transcript` | Only complete transcripts will be considered for SV annotation. Default: True. |
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
  | 14 | affinity | Binding affinity (nM) provided by NetMHCpan |
  | 15 | BA_rank | BA rank of the binding provided by NetMHCpan |
  | 16 | EL_rank | EL rank of the binding provided by NetMHCpan. From NetMHCpan4.0, EL rank is the most recommended feature for filtering neoantigens.|

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
