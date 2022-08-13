import re
from collections import OrderedDict
from neosv.sv_class import StructuralVariant


def get_coordinate(alt):
    """
    :param alt: the ALT field of VCF format
    :return: chrom and position of this breakpoint
    """
    coord_pattern = re.compile(r'(\]|\[)(.+?):(\d+)(\]|\[)')
    coord = coord_pattern.findall(alt)[0]
    chrom, pos = coord[1].replace('chr', ''), int(coord[2])
    return chrom, pos


def get_insertion(alt):
    """
    :param alt: the ALT field of VCF format
    :return: inserted sequence of this SV
    """
    nt_pattern = re.compile(r'[ATCG]+')
    right_pattern = re.compile(r'([ATCG]+)(\]|\[)')
    nt = nt_pattern.findall(alt)[0]
    if len(nt) == 1:
        return ''
    else:
        if right_pattern.search(alt):
            return nt[1:]
        else:
            return nt[:-1]


def sv_pattern_infer(svvcf):
    """
    :param svvcf: a SV stored as a VariantCallingFormat object
    :return: a StructuralVariant object
    """
    chrom1, pos1 = svvcf.chrom.replace('chr', ''), int(svvcf.pos)
    chrom2, pos2 = get_coordinate(svvcf.alt)
    insertion = get_insertion(svvcf.alt)
    # t[p[
    pattern1 = re.compile(r'[ATCG]+\[.+?\[')
    # t]p]
    pattern2 = re.compile(r'[ATCG]+\].+?\]')
    # ]p]t
    pattern3 = re.compile(r'\].+?\][ATCG]+')
    # [p[t
    pattern4 = re.compile(r'\[.+?\[[ATCG]+')

    if pattern1.search(svvcf.alt):
        pattern = 1
    elif pattern2.search(svvcf.alt):
        pattern = 2
    elif pattern3.search(svvcf.alt):
        pattern = 3
    else:
        pattern = 4
    return StructuralVariant(chrom1, pos1, chrom2, pos2, insertion, pattern)


def remove_duplicate(sv_list):
    """
    :param sv_list: a list of StructuralVariant objects
    :return: deduplicated sv_list, maintaining the initial order
    """
    ordered_dict = OrderedDict.fromkeys(sv_list)
    sv_list_rmdup = list(ordered_dict.keys())
    print("{0} duplicated SVs were removed".format(len(sv_list)-len(sv_list_rmdup)))
    return sv_list_rmdup
