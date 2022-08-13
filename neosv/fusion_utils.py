from neosv.fusion_class import CDS, CDSCollection, SVFusion
from neosv.transcript_utils import get_cds_range, get_noncds_range, is_overlap, get_transcript


def truncate_cds(transcript, direction, pos):
    cds_ranges = get_cds_range(transcript)
    ncds_ranges = get_noncds_range(transcript)
    cds_overlap, cds_index = is_overlap(pos, cds_ranges)
    ncds_overlap, ncds_index = is_overlap(pos, ncds_ranges)
    if transcript.strand == '+' and direction == '5':
        if cds_overlap:
            cds_list = [CDS(cds[0],
                            cds[1],
                            True,
                            False,
                            False
                            ) for cds in cds_ranges[: cds_index+1]]
            cds_list[0].startcodon = True
            cds_list[-1].end = pos
            cds_list[-1].intact = False
        else:
            cds_list = [CDS(cds[0],
                            cds[1],
                            True,
                            False,
                            False
                            ) for cds in cds_ranges[: ncds_index]]
            if cds_list:
                cds_list[0].startcodon = True
            if ncds_index == len(cds_ranges):
                cds_list[-1].stopcodon = True
    elif transcript.strand == '+' and direction == '3':
        if cds_overlap:
            cds_list = [CDS(cds[0],
                            cds[1],
                            True,
                            False,
                            False
                            ) for cds in cds_ranges[cds_index: ]]
            cds_list[0].start = pos
            cds_list[0].intact = False
            cds_list[-1].stopcodon = True
        else:
            cds_list = [CDS(cds[0],
                            cds[1],
                            True,
                            False,
                            False
                            ) for cds in cds_ranges[ncds_index: ]]
            if cds_list:
                cds_list[-1].stopcodon = True
            if ncds_index == 0:
                cds_list[0].startcodon = True
    elif transcript.strand == '-' and direction == '5':
        if cds_overlap:
            cds_list = [CDS(cds[0],
                            cds[1],
                            True,
                            False,
                            False
                            ) for cds in cds_ranges[: cds_index+1]]
            cds_list[0].startcodon = True
            cds_list[-1].start = pos
            cds_list[-1].intact = False
        else:
            cds_list = [CDS(cds[0],
                            cds[1],
                            True,
                            False,
                            False
                            ) for cds in cds_ranges[: ncds_index]]
            if cds_list:
                cds_list[0].startcodon = True
            if ncds_index == len(cds_ranges):
                cds_list[-1].stopcodon = True
    else:
        if cds_overlap:
            cds_list = [CDS(cds[0],
                            cds[1],
                            True,
                            False,
                            False
                            ) for cds in cds_ranges[cds_index:]]
            cds_list[0].start = pos
            cds_list[0].intact = False
            cds_list[-1].stopcodon = True
        else:
            cds_list = [CDS(cds[0],
                            cds[1],
                            True,
                            False,
                            False
                            ) for cds in cds_ranges[ncds_index:]]
            if cds_list:
                cds_list[-1].stopcodon = True
            if ncds_index == 0:
                cds_list[0].startcodon = True
    return CDSCollection(transcript, cds_list, direction)


def sv_to_svfusion(sv, ensembl):
    transcript_1 = get_transcript(sv.chrom1, sv.pos1, ensembl, complete=True)
    transcript_2 = get_transcript(sv.chrom2, sv.pos2, ensembl, complete=True)
    if transcript_1 and transcript_2:
        if transcript_1.strand == '+' and transcript_2.strand == '+':
            if sv.pattern == 1:
                cdscollection_1 = truncate_cds(transcript_1, '5', sv.pos1)
                cdscollection_2 = truncate_cds(transcript_2, '3', sv.pos2)
            elif sv.pattern == 2:
                cdscollection_1 = None
                cdscollection_2 = None
            elif sv.pattern == 3:
                cdscollection_1 = truncate_cds(transcript_1, '3', sv.pos1)
                cdscollection_2 = truncate_cds(transcript_2, '5', sv.pos2)
            else:
                cdscollection_1 = None
                cdscollection_2 = None
        elif transcript_1.strand == '+' and transcript_2.strand == '-':
            if sv.pattern == 1:
                cdscollection_1 = None
                cdscollection_2 = None
            elif sv.pattern == 2:
                cdscollection_1 = truncate_cds(transcript_1, '5', sv.pos1)
                cdscollection_2 = truncate_cds(transcript_2, '3', sv.pos2)
            elif sv.pattern == 3:
                cdscollection_1 = None
                cdscollection_2 = None
            else:
                cdscollection_1 = truncate_cds(transcript_1, '3', sv.pos1)
                cdscollection_2 = truncate_cds(transcript_2, '5', sv.pos2)
        elif transcript_1.strand == '-' and transcript_2.strand == '+':
            if sv.pattern == 1:
                cdscollection_1 = None
                cdscollection_2 = None
            elif sv.pattern == 2:
                cdscollection_1 = truncate_cds(transcript_1, '3', sv.pos1)
                cdscollection_2 = truncate_cds(transcript_2, '5', sv.pos2)
            elif sv.pattern == 3:
                cdscollection_1 = None
                cdscollection_2 = None
            else:
                cdscollection_1 = truncate_cds(transcript_1, '5', sv.pos1)
                cdscollection_2 = truncate_cds(transcript_2, '3', sv.pos2)
        else:
            if sv.pattern == 1:
                cdscollection_1 = truncate_cds(transcript_1, '3', sv.pos1)
                cdscollection_2 = truncate_cds(transcript_2, '5', sv.pos2)
            elif sv.pattern == 2:
                cdscollection_1 = None
                cdscollection_2 = None
            elif sv.pattern == 3:
                cdscollection_1 = truncate_cds(transcript_1, '5', sv.pos1)
                cdscollection_2 = truncate_cds(transcript_2, '3', sv.pos2)
            else:
                cdscollection_1 = None
                cdscollection_2 = None
    else:
        cdscollection_1 = None
        cdscollection_2 = None
    return SVFusion(sv, cdscollection_1, cdscollection_2)


def svfusions_to_dict(svfusions):
    """
    :param svfusions: a list of svfusion class
    :return: a dict whose key is svfusion, value is neoepitopes
    """
    svfusion_dic = {}
    for svfusion in svfusions:
        svfusion_dic[svfusion] = svfusion.neoepitopes
    return svfusion_dic
