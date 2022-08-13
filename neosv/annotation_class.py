class EXON(object):
    """
    Class for the genomic coordinates of an
    exon, and whether it's an intact exon.
    NOTE: start is always smaller than end
    """
    def __init__(self, start, end, intact):
        self.start = int(start)
        self.end = int(end)
        self.intact = intact

    def __str__(self):
        return "%s(start = %d, end = %d, intact = %r)" % (
            self.__class__.__name__,
            self.start,
            self.end,
            self.intact
        )

    def __repr__(self):
        return "%s(%d, %d, %r)" % (
            self.__class__.__name__,
            self.start,
            self.end,
            self.intact
        )


class EXONCollection(object):
    """
    Class for a collection of EXONs
    transcript: Transcript class of pyensembl
    exon_list: a list of EXONs
    exon_idx: index of the exon in which the transcript end, if end in intron, -1
    intron_idx: index of the exon in which the transcript end, if end in exon, -1
    part: which part of the transcript is retained, 3 or 5 prime
    """
    def __init__(self, transcript, exon_list, exon_idx, intron_idx, part):
        self.transcript = transcript
        self.exon_list = exon_list
        self.exon_idx = exon_idx
        self.intron_idx = intron_idx
        self.part = part

    def __str__(self):
        return "%s(transcript_id = %s, transcript_name = %s, exon_idx = %d, intron_idx = %d, part= %s)" % (
            self.__class__.__name__,
            self.transcript.transcript_id,
            self.transcript.transcript_name,
            self.exon_idx,
            self.intron_idx,
            self.part
        )

    def __repr__(self):
        return "%s(%s, %s, %d, %d, %s)" % (
            self.__class__.__name__,
            self.transcript.transcript_id,
            self.transcript.transcript_name,
            self.exon_idx,
            self.intron_idx,
            self.part
        )

    @property
    def transcript_id(self):
        """
        :return: transcript id of this collection of exons
        """
        return self.transcript.transcript_id

    @property
    def gene_name(self):
        """
        :return: gene name of this collection of exons
        """
        return self.transcript.gene.gene_name

    @property
    def transcript_name(self):
        """
        :return: transcript name of this collection of exons
        """
        return self.transcript.transcript_name

    @property
    def exon_num(self):
        """
        :return: number of exons in this collection
        """
        return len(self.transcript.exons)

    @property
    def breakloc(self):
        """
        :return: which element this collection ends, exon or intron
        """
        if self.exon_idx == -1:
            return 'Intron'
        else:
            return 'Exon'

    @property
    def _transcript_retain(self):
        """
        :return: a list of the elements retained in this collection,
                 the element containing the breakend is in lower case
        """
        if self.part == '5' and self.breakloc == 'Intron':
            transcript_retain = ['E{}-I{}'.format(i+1, i+1) for i in range(0, self.intron_idx+1)]
            transcript_retain[-1] = transcript_retain[-1].replace('I', 'i')
        elif self.part == '5' and self.breakloc == 'Exon':
            transcript_retain = ['E{}-I{}'.format(i+1, i+1) for i in range(0, self.exon_idx)] + \
                                ['e{}'.format(self.exon_idx+1)]
        elif self.part == '3' and self.breakloc == 'Intron':
            transcript_retain = ['I{}-E{}'.format(i+1, i+2) for i in range(self.intron_idx, self.exon_num-1)]
            transcript_retain[0] = transcript_retain[0].replace('I', 'i')
        else:
            transcript_retain = ['e{}'.format(self.exon_idx+1)] + \
                                ['I{}-E{}'.format(i+1, i+2) for i in range(self.exon_idx, self.exon_num-1)]
        return transcript_retain

    def annot_l(self):
        """
        :return: the long format of transcript
        """
        annot = '-'.join(self._transcript_retain)
        return annot

    def annot_s(self):
        """
        :return: the short format of transcript
        """
        annot = '-'.join(self._transcript_retain)
        annot_split = annot.split('-')
        if len(annot_split) <= 2:
            return '-'.join(annot_split)
        else:
            return '-'.join([annot_split[0], annot_split[-1]])


class SVEffect(object):
    """
    Class for annotating the effect of a SV
    sv: StructuralVariant class
    ec_1 and ec_2: EXONCollection class
    """
    def __init__(self, sv, ec_1, ec_2):
        self.sv = sv
        self.ec_1 = ec_1
        self.ec_2 = ec_2

    def __str__(self):
        return "%s(sv = %r, ec_1 = %r, ec_2 = %r)" % (
            self.__class__.__name__,
            self.sv,
            self.ec_1,
            self.ec_2
        )

    def __repr__(self):
        return "%s(%r, %r, %r)" % (
            self.__class__.__name__,
            self.sv,
            self.ec_1,
            self.ec_2
        )

    @property
    def fusion(self):
        if self.ec_1 and self.ec_2:
            if self.ec_1.transcript.strand == self.ec_2.transcript.strand:
                if self.sv.pattern in [1, 3]:
                    return True
                else:
                    return False
            else:
                if self.sv.pattern in [1, 3]:
                    return False
                else:
                    return True
        else:
            return False

    def output(self, annoformat='short'):
        tran_annot_1, tran_annot_2 = 'None', 'None'
        tran_id_1, tran_id_2 = 'None', 'None'
        gene_name_1, gene_name_2 = 'None', 'None'
        loc_1, loc_2 = 'Intergenic', 'Intergenic'
        strand_1, strand_2 = 'None', 'None'
        fusion = 'No-fusion'
        if annoformat == 'long':
            if self.ec_1:
                tran_id_1 = self.ec_1.transcript_id
                tran_annot_1 = self.ec_1.annot_l()
                gene_name_1 = self.ec_1.gene_name
                loc_1 = self.ec_1.breakloc
                strand_1 = self.ec_1.transcript.strand
            if self.ec_2:
                tran_id_2 = self.ec_2.transcript_id
                tran_annot_2 = self.ec_2.annot_l()
                gene_name_2 = self.ec_2.gene_name
                loc_2 = self.ec_2.breakloc
                strand_2 = self.ec_2.transcript.strand
        else:
            if self.ec_1:
                tran_id_1 = self.ec_1.transcript_id
                tran_annot_1 = self.ec_1.annot_s()
                gene_name_1 = self.ec_1.gene_name
                loc_1 = self.ec_1.breakloc
                strand_1 = self.ec_1.transcript.strand
            if self.ec_2:
                tran_id_2 = self.ec_2.transcript_id
                tran_annot_2 = self.ec_2.annot_s()
                gene_name_2 = self.ec_2.gene_name
                loc_2 = self.ec_2.breakloc
                strand_2 = self.ec_2.transcript.strand

        chrom_1, chrom_2 = self.sv.chrom1, self.sv.chrom2
        pos_1, pos_2 = self.sv.pos1, self.sv.pos2
        if self.fusion:
            fusion = 'Fusion'
        return [chrom_1, str(pos_1), loc_1, gene_name_1, tran_id_1, strand_1, tran_annot_1,
                chrom_2, str(pos_2), loc_2, gene_name_2, tran_id_2, strand_2, tran_annot_2,
                str(self.sv.pattern), self.sv.svtype, fusion]
