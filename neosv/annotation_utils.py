from neosv.annotation_class import EXON, EXONCollection, SVEffect
from neosv.transcript_utils import get_transcript, get_exon_range, get_intron_range, is_overlap


def truncate_exon(transcript, direction, pos):
    """
    :param transcript: a Transcript class of pyensembl
    :param direction: 5(upstream) or 3(downstream) of the breakpoint
    :param pos: position of the breakpoint
    :return: an EXONCollection class
    """
    exon_ranges = get_exon_range(transcript)
    intron_ranges = get_intron_range(transcript)
    exon_overlap, exon_index = is_overlap(pos, exon_ranges)
    intron_overlap, intron_index = is_overlap(pos, intron_ranges)
    if transcript.strand == '+' and direction == '5':
        if exon_overlap:
            exon_list = [EXON(exon[0],
                              exon[1],
                              True
                              ) for exon in exon_ranges[:exon_index+1]]
            exon_list[-1].end = pos
            exon_list[-1].intact = False
        else:
            exon_list = [EXON(exon[0],
                              exon[1],
                              True
                              ) for exon in exon_ranges[:intron_index+1]]
    elif transcript.strand == '+' and direction == '3':
        if exon_overlap:
            exon_list = [EXON(exon[0],
                              exon[1],
                              True
                              ) for exon in exon_ranges[exon_index:]]
            exon_list[0].start = pos
            exon_list[0].intact = False
        else:
            exon_list = [EXON(exon[0],
                              exon[1],
                              True
                              ) for exon in exon_ranges[intron_index+1:]]
    elif transcript.strand == '-' and direction == '5':
        if exon_overlap:
            exon_list = [EXON(exon[0],
                              exon[1],
                              True
                              ) for exon in exon_ranges[: exon_index+1]]
            exon_list[-1].start = pos
            exon_list[-1].intact = False
        else:
            exon_list = [EXON(exon[0],
                              exon[1],
                              True
                              ) for exon in exon_ranges[:intron_index+1]]
    else:
        if exon_overlap:
            exon_list = [EXON(exon[0],
                              exon[1],
                              True
                              ) for exon in exon_ranges[exon_index:]]
            exon_list[0].end = pos
            exon_list[0].intact = False
        else:
            exon_list = [EXON(exon[0],
                              exon[1],
                              True
                              ) for exon in exon_ranges[intron_index+1:]]
    return EXONCollection(transcript, exon_list, exon_index, intron_index, direction)


def sv_to_sveffect(sv, ensembl, complete=False):
    """
    :param sv: a StructuralVariant class
    :param ensembl: a Genome class in pyensembl
    :param complete: whether to only use complete transcripts
    :return: SVEffect class
    """
    transcript_1 = get_transcript(sv.chrom1, sv.pos1, ensembl, complete)
    transcript_2 = get_transcript(sv.chrom2, sv.pos2, ensembl, complete)
    if transcript_1 and transcript_2:
        if transcript_1.strand == '+' and transcript_2.strand == '+':
            if sv.pattern == 1:
                exoncollection_1 = truncate_exon(transcript_1, '5', sv.pos1)
                exoncollection_2 = truncate_exon(transcript_2, '3', sv.pos2)
                return [SVEffect(sv, exoncollection_1, exoncollection_2)]
            elif sv.pattern == 2:
                exoncollection_1 = truncate_exon(transcript_1, '5', sv.pos1)
                exoncollection_2 = truncate_exon(transcript_2, '5', sv.pos2)
                return [SVEffect(sv, exoncollection_1, None),
                        SVEffect(sv, None, exoncollection_2)]
            elif sv.pattern == 3:
                exoncollection_1 = truncate_exon(transcript_1, '3', sv.pos1)
                exoncollection_2 = truncate_exon(transcript_2, '5', sv.pos2)
                return [SVEffect(sv, exoncollection_1, exoncollection_2)]
            else:
                exoncollection_1 = truncate_exon(transcript_1, '3', sv.pos1)
                exoncollection_2 = truncate_exon(transcript_2, '3', sv.pos2)
                return [SVEffect(sv, exoncollection_1, None),
                        SVEffect(sv, None, exoncollection_2)]
        elif transcript_1.strand == '+' and transcript_2.strand == '-':
            if sv.pattern == 1:
                exoncollection_1 = truncate_exon(transcript_1, '5', sv.pos1)
                exoncollection_2 = truncate_exon(transcript_2, '5', sv.pos2)
                return [SVEffect(sv, exoncollection_1, None),
                        SVEffect(sv, None, exoncollection_2)]
            elif sv.pattern == 2:
                exoncollection_1 = truncate_exon(transcript_1, '5', sv.pos1)
                exoncollection_2 = truncate_exon(transcript_2, '3', sv.pos2)
                return [SVEffect(sv, exoncollection_1, exoncollection_2)]
            elif sv.pattern == 3:
                exoncollection_1 = truncate_exon(transcript_1, '3', sv.pos1)
                exoncollection_2 = truncate_exon(transcript_2, '3', sv.pos2)
                return [SVEffect(sv, exoncollection_1, None),
                        SVEffect(sv, None, exoncollection_2)]
            else:
                exoncollection_1 = truncate_exon(transcript_1, '3', sv.pos1)
                exoncollection_2 = truncate_exon(transcript_2, '5', sv.pos2)
                return [SVEffect(sv, exoncollection_1, exoncollection_2)]
        elif transcript_1.strand == '-' and transcript_2.strand == '+':
            if sv.pattern == 1:
                exoncollection_1 = truncate_exon(transcript_1, '3', sv.pos1)
                exoncollection_2 = truncate_exon(transcript_2, '3', sv.pos2)
                return [SVEffect(sv, exoncollection_1, None),
                        SVEffect(sv, None, exoncollection_2)]
            elif sv.pattern == 2:
                exoncollection_1 = truncate_exon(transcript_1, '3', sv.pos1)
                exoncollection_2 = truncate_exon(transcript_2, '5', sv.pos2)
                return [SVEffect(sv, exoncollection_1, exoncollection_2)]
            elif sv.pattern == 3:
                exoncollection_1 = truncate_exon(transcript_1, '5', sv.pos1)
                exoncollection_2 = truncate_exon(transcript_2, '5', sv.pos2)
                return [SVEffect(sv, exoncollection_1, None),
                        SVEffect(sv, None, exoncollection_2)]
            else:
                exoncollection_1 = truncate_exon(transcript_1, '5', sv.pos1)
                exoncollection_2 = truncate_exon(transcript_2, '3', sv.pos2)
                return [SVEffect(sv, exoncollection_1, exoncollection_2)]
        else:
            if sv.pattern == 1:
                exoncollection_1 = truncate_exon(transcript_1, '3', sv.pos1)
                exoncollection_2 = truncate_exon(transcript_2, '5', sv.pos2)
                return [SVEffect(sv, exoncollection_1, exoncollection_2)]
            elif sv.pattern == 2:
                exoncollection_1 = truncate_exon(transcript_1, '3', sv.pos1)
                exoncollection_2 = truncate_exon(transcript_2, '3', sv.pos2)
                return [SVEffect(sv, exoncollection_1, None),
                        SVEffect(sv, None, exoncollection_2)]
            elif sv.pattern == 3:
                exoncollection_1 = truncate_exon(transcript_1, '5', sv.pos1)
                exoncollection_2 = truncate_exon(transcript_2, '3', sv.pos2)
                return [SVEffect(sv, exoncollection_1, exoncollection_2)]
            else:
                exoncollection_1 = truncate_exon(transcript_1, '5', sv.pos1)
                exoncollection_2 = truncate_exon(transcript_2, '5', sv.pos2)
                return [SVEffect(sv, exoncollection_1, None),
                        SVEffect(sv, None, exoncollection_2)]
    elif transcript_1 and not transcript_2:
        if transcript_1.strand == '+':
            if sv.pattern in [1, 2]:
                exoncollection = truncate_exon(transcript_1, '5', sv.pos1)
                return [SVEffect(sv, exoncollection, None)]
            else:
                exoncollection = truncate_exon(transcript_1, '3', sv.pos1)
                return [SVEffect(sv, exoncollection, None)]
        else:
            if sv.pattern in [1, 2]:
                exoncollection = truncate_exon(transcript_1, '3', sv.pos1)
                return [SVEffect(sv, exoncollection, None)]
            else:
                exoncollection = truncate_exon(transcript_1, '5', sv.pos1)
                return [SVEffect(sv, exoncollection, None)]
    elif not transcript_1 and transcript_2:
        if transcript_2.strand == '+':
            if sv.pattern in [1, 4]:
                exoncollection = truncate_exon(transcript_2, '3', sv.pos2)
                return [SVEffect(sv, None, exoncollection)]
            else:
                exoncollection = truncate_exon(transcript_2, '5', sv.pos2)
                return [SVEffect(sv, None, exoncollection)]
        else:
            if sv.pattern in [1, 4]:
                exoncollection = truncate_exon(transcript_2, '5', sv.pos2)
                return [SVEffect(sv, None, exoncollection)]
            else:
                exoncollection = truncate_exon(transcript_2, '3', sv.pos2)
                return [SVEffect(sv, None, exoncollection)]
    else:
        return [SVEffect(sv, None, None)]
