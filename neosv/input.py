import os
import re
from warnings import warn
from pyensembl.genome import Genome
from pyensembl import EnsemblRelease


class VariantCallingFormat(object):
    """
    Class for storing SV information in VCF format,
    all components are in string format
    """

    def __init__(self, chrom, pos, ref, alt):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt

    def __str__(self):
        return "%s(chrom = %s, pos = %s, ref = %s, alt = %s)" % (
            self.__class__.__name__,
            self.chrom,
            self.pos,
            self.ref,
            self.alt
        )

    def __repr__(self):
        return "%s(%s, %s, %s, %s)" % (
            self.__class__.__name__,
            self.chrom,
            self.pos,
            self.ref,
            self.alt
        )


class BedpeFormat(object):
    """
    Class for storing SV information in BEDPE format,
    all components are in string format
    """

    def __init__(self, chrom1, pos1, strand1, chrom2, pos2, strand2):
        self.chrom1 = chrom1
        self.pos1 = pos1
        self.strand1 = strand1
        self.chrom2 = chrom2
        self.pos2 = pos2
        self.strand2 = strand2

    def __str__(self):
        return "%s(chrom1 = %s, pos1 = %s, strand1 = %s, chrom2 = %s, pos2 = %s, strand2 = %s)" % (
            self.__class__.__name__,
            self.chrom1,
        	self.pos1,
        	self.strand1,
        	self.chrom2,
        	self.pos2,
        	self.strand2
        )

    def __repr__(self):
        return "%s(%s, %s, %s, %s, %s, %s)" % (
            self.__class__.__name__,
            self.chrom1,
        	self.pos1,
        	self.strand1,
        	self.chrom2,
        	self.pos2,
        	self.strand2
        )


def vcf_load(filepath):
    """
    :param filepath: the absolute path of a VCF file
    :return: a list of VariantCallingFormat objects
    """
    svvcf_list = []
    filename = os.path.basename(filepath)
    line_num = 0
    print("Loading SVs from {0}.".format(filename))
    with open(filepath, 'r') as f:
        for line in f:
            line_num += 1
            if line.startswith('#'):
                continue
            else:
                tmpline = line.rstrip().split("\t")
                chrom = tmpline[0].replace('chr', '')
                pos = tmpline[1]
                ref = tmpline[3]
                alt = tmpline[4].replace('chr', '')
                if vcf_alt_format_check(alt):
                    svvcf = VariantCallingFormat(chrom, pos, ref, alt)
                    svvcf_list.append(svvcf)
                else:
                    warn("File:\t{0}\tLine:\t{1}\nMalformed SV:\t{2}\t{3}\t{4}\t{5}".format(
                        filename, line_num, chrom, pos, ref, alt
                    ))
    if not svvcf_list:
        warn("No SV was detected in {0}.".format(filename))
    else:
        print("{0} SVs were detected.".format(len(svvcf_list)))
    return svvcf_list


def vcf_alt_format_check(alt):
    """
    :param alt: the ALT field of a VCF file
    :return: bool value indicate whether it is a legal ALT field
    """
    legal_pattern_right = re.compile(r'([ATCGatcg]+)(\]|\[)(.+?):(\d+)(\]|\[)$')
    legal_pattern_left = re.compile(r'(\]|\[)(.+?):(\d+)(\]|\[)([ATCGatcg]+)$')
    isright = legal_pattern_right.match(alt)
    isleft = legal_pattern_left.match(alt)
    return isright or isleft


def bedpe_load(filepath):
    """
    :param filepath: the absolute path of a BEDPE file
    :return: a list of BEDPE objects
    """
    bedpe_list = []
    filename = os.path.basename(filepath)
    line_num = 0
    print("Loading SVs from {0}.".format(filename))
    with open(filepath, 'r') as f:
        header = next(f)
        header = header.rstrip().split('\t')
        for line in f:
            line_num += 1
            tmpline = line.rstrip().split("\t")
            chrom1 = tmpline[header.index('chrom1')].replace('chr', '')
            pos1 = tmpline[header.index('start1')]
            chrom2 = tmpline[header.index('chrom2')].replace('chr', '')
            pos2 = tmpline[header.index('start2')]
            strand1 = tmpline[header.index('strand1')]
            strand2 = tmpline[header.index('strand2')]
            bedpe = BedpeFormat(chrom1, pos1, strand1, chrom2, pos2, strand2)
            bedpe_list.append(bedpe)

    if not bedpe_list:
        warn("No SV was detected in {0}.".format(filename))
    else:
        print("{0} SVs were detected.".format(len(bedpe_list)))
    return bedpe_list


def hla_load(filepath):
    """
    :param filepath: the absolute path of a HLA typing file
    :return: a list of HLA alleles joined by ,
    """
    hla_alleles = []
    filename = os.path.basename(filepath)
    with open(filepath, 'r') as f:
        for line in f:
            hla_allele = line.rstrip()
            hla_alleles.append(hla_allele.replace('*', ''))
    return ','.join(hla_alleles)



def ensembl_load(release, gtf_file, cdna_file, cache_dir):
    """
    :param release: the release number in EMSEMBL, could be custom
    :param gtf_file: the path of gtf file if release == custom
    :param cdna_file: the path of cdna file if release == custom
    :param cache_dir: directory for pyensembl downloading
    :return: a Genome class in pyensembl
    """
    if cache_dir:
        os.environ['PYENSEMBL_CACHE_DIR'] = cache_dir
    if release != 'custom':
        ensembl = EnsemblRelease(int(release))
        ensembl.download()
        ensembl.index()
    else:
        ensembl = Genome(gtf_path_or_url=gtf_file,
                         transcript_fasta_paths_or_urls=cdna_file,
                         reference_name='User-defined',
                         annotation_name='User-defined')
        ensembl.index()
    return ensembl


def get_window_range(window):
    """
    :param window: the args.window given by user
    :return: a list of window sizes
    """
    start = int(window.split('-')[0])
    end = int(window.split('-')[1])
    return list(range(start, end+1))