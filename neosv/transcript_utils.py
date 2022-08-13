def get_longest_transcript(transcripts):
    """
    :param transcripts: a list of Transcript(pyensembl) instances
    :return: the longest transcript
    """
    transcripts = sorted(transcripts, key=lambda t: t.end-t.start)
    transcript = transcripts[-1]
    return transcript


def get_transcript(chrom, pos, ensembl, complete=True):
    """
    :param chrom: chromosome with no chr of breakpoint
    :param pos: position of breakpoint
    :param ensembl: Genome instance in pyensembl
    :param complete: only consider complete transcripts
    :return: firstly return the longest complete transcript,
             if there is no complete transcript and
             complete = False, return the longest transcript
    """
    transcripts = ensembl.transcripts_at_locus(contig=str(chrom), position=int(pos))
    transcripts_comp = [transcript for transcript in transcripts if transcript.complete]
    if transcripts_comp:
        return get_longest_transcript(transcripts_comp)
    else:
        if complete:
            return None
        else:
            if transcripts:
                return get_longest_transcript(transcripts)
            else:
                return None


def get_cds_range(transcript):
    """
    :param transcript: transcript instance in pyensembl
    :return: cds intervals of this transcript
             from 5' to 3', cds 1, cds 2, ...
             [start, end], start < end
    """
    cds_ranges = transcript.coding_sequence_position_ranges
    return cds_ranges


def get_noncds_range(transcript):
    """
    :param transcript: transcript instance in pyensembl
    :return: noncds intervals of this transcript
             from 5' to 3', noncds 1, noncds 2, ...
             [start, end], start < end
             if there is no interval, for example, no utr,
             we will add a interval (0,0)
    """
    cds_ranges = transcript.coding_sequence_position_ranges
    ncds_ranges = []
    if transcript.strand == '+':
        utr5_start = transcript.start
        utr5_end = cds_ranges[0][0] - 1
        if utr5_start > utr5_end:
            ncds_ranges.append((0, 0))
        else:
            ncds_ranges.append((utr5_start, utr5_end))
        for i in range(1, len(cds_ranges)):
            cds_start = cds_ranges[i-1][1] + 1
            cds_end = cds_ranges[i][0] - 1
            ncds_ranges.append((cds_start, cds_end))
        utr3_start = cds_ranges[-1][1] + 1
        utr3_end = transcript.end
        if utr3_start > utr3_end:
            ncds_ranges.append((0, 0))
        else:
            ncds_ranges.append((utr3_start, utr3_end))
    else:
        utr5_start = cds_ranges[0][1] + 1
        utr5_end = transcript.end
        if utr5_start > utr5_end:
            ncds_ranges.append((0, 0))
        else:
            ncds_ranges.append((utr5_start, utr5_end))
        for i in range(1, len(cds_ranges)):
            cds_start = cds_ranges[i][1] + 1
            cds_end = cds_ranges[i-1][0] - 1
            ncds_ranges.append((cds_start, cds_end))
        utr3_start = transcript.start
        utr3_end = cds_ranges[-1][0] - 1
        if utr3_start > utr3_end:
            ncds_ranges.append((0, 0))
        else:
            ncds_ranges.append((utr3_start, utr3_end))
    return ncds_ranges


def get_exon_range(transcript):
    """
    :param transcript: transcript instance in pyensembl
    :return: exon intervals of this transcript
             from 5' to 3', exon 1, exon 2, ...
             [start, end], start < end
    """
    exon_ranges = []
    for exon in transcript.exons:
        exon_ranges.append((exon.start, exon.end))
    return exon_ranges


def get_intron_range(transcript):
    """
    :param transcript: transcript instance in pyensembl
    :return: non-exon intervals of this transcript
             from 5' to 3', intron 1, intron 2, ...
             [start, end], start < end
    """
    exon_ranges = get_exon_range(transcript)
    intron_ranges = []
    if transcript.strand == '+':
        for i in range(1, len(exon_ranges)):
            intron_start = exon_ranges[i-1][1] + 1
            intron_end = exon_ranges[i][0] - 1
            intron_ranges.append((intron_start, intron_end))
    else:
        for i in range(1, len(exon_ranges)):
            intron_start = exon_ranges[i][1] + 1
            intron_end = exon_ranges[i-1][0] - 1
            intron_ranges.append((intron_start, intron_end))
    if intron_ranges:
        return intron_ranges
    else:
        return [(0, 0)]


def is_overlap(pos, intervals):
    """
    :param pos: a genomic position
    :param intervals: a list of genomic intervals
    :return: whether the pos overlaps with one of intervals
             and the index of overlapped interval
    """
    for i in range(len(intervals)):
        if intervals[i][0] <= pos <= intervals[i][1]:
            return True, i
    return False, -1
