import argparse
import os
import sys


def create_arg_parser():
    parser = argparse.ArgumentParser(prog="neosv")
    parser.add_argument('-vf', '--vcf-file', dest='vcffile', metavar='VCF_FILE', required=True,
                        help='Structural variants in VCF format..')
    parser.add_argument('-hf', '--hla-file', dest='hlafile', metavar='HLA_FILE', default=None,
                        help='HLA alleles (resolution: 4 digit, e.g. HLA-A*02:01), with one allele per line.')
    parser.add_argument('-np', '--netmhc-path', dest='netmhc', metavar='NETMHC_PATH', default=None,
                        help='Absolute path to the netMHCpan execution file, please skip this parameter if netMHCpan has been added to your PATH.')
    parser.add_argument('-o', '--out', dest='outdir', metavar='OUTDIR', required=True,
                        help="Folder for all result files. A new folder will be created if it does not exist.")
    parser.add_argument('-p', '--prefix', dest='prefix', metavar='PREFIX', default='sample',
                        help="This prefix will be added to all output files.")
    parser.add_argument('-r', '--release', dest='release', metavar='RELEASE', default='75',
                        help='Which reference (ENSEMBL release) you want to use. Ensembl releases that'
                             'correspond to hg18/NCBI36, hg19/GRCh37, hg38/GRCh38 are 54, 75, 95.'
                             'If your data are from other species(custom), please download the gtf '
                             'file and the cdna file from ENSEMBL website ftp://ftp.ensembl.org/pub'
                             ' and specify them using --gtf-file and --cdna-file.')
    parser.add_argument('-gf', '--gtf-file', dest='gtffile', metavar='GTF_FILE', default=None,
                        help='GTF file for the reference.')
    parser.add_argument('-cf', '--cdna-file', dest='cdnafile', metavar='CDNA_FILE', default=None,
                        help='cDNA file for the reference.')
    parser.add_argument('-pd', '--pyensembl-cache-dir', dest='cachedir', metavar='PYEMSEMBL_CACHE_DIR', default=None,
                        help='Directory for Pyensembl cache files. If not specified, the platform-specific cache folder will be used.')
    parser.add_argument('-l', '--epitope-lengths', dest='window', metavar='EPITOPE_LENGTHS', default='8-11',
                        help='Lengths of neoepitopes to predict MHC binding. Default: 8-11.')
    parser.add_argument('-ic', '--ic50-cutoff', dest='aff_cutoff', metavar='IC50_CUTOFF', type=float, default=500,
                        help='Filter neoepitopes with IC50 (nM) above this value. Default: 500.')
    parser.add_argument('-rc', '--ranking-cutoff', dest='rank_cutoff', metavar='RANKING_CUTOFF', type=float,
                        default='2.0', help='Filter neoepitopes with rank above this value. Default: 2.')
    parser.add_argument('-ct', '--complete-transcript', dest='complete', metavar='COMPLETE_TRANSCRIPT', type=bool,
                        default=False, help='Only complete transcripts will be considered for SV annotation.')
    parser.add_argument('--anno-only', dest='anno', action='store_true', default=False,
                        help='Only annotate SV without predicting neoantigens.If this argument is added,'
                        '--hla-file is not required, and you will only get the annotation result.')

    args = parser.parse_args()

    if args.release != 'custom':
        valid_range = [str(i) for i in range(54, 107)]
        if args.release not in valid_range:
            sys.exit('Release number must be between 54 and 107.')
    elif args.release == 'custom':
        if not args.gtffile or not args.cdnafile:
            sys.exit('A custom release is used, but no gtffile and cdnafile specified.')
    else:
        sys.exit('--release must be a number between 54-107 or custom.')

    if not args.anno and not args.hlafile:
        sys.exit('No HLA file specified. The annotation-only mode could be used '
                 '(--annotation-only) if you do not have HLA information.')
    if not args.anno and not args.netmhc:
        sys.exit('No netMHC exec file specified. he annotation-only mode could be used '
                 '(--annotation-only) if you do not want to predict neoantigen.')
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)

    return args
