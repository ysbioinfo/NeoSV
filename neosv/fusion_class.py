from Bio.Seq import Seq


class CDS(object):
    def __init__(self, start, end, intact, startcodon, stopcodon):
        self.start = int(start)
        self.end = int(end)
        self.intact = intact
        self.startcodon = startcodon
        self.stopcodon = stopcodon

    def __str__(self):
        return "%s(start = %d, end = %d, intact = %r, startcodon = %r, stopcodon = %r)" % (
            self.__class__.__name__,
            self.start,
            self.end,
            self.intact,
            self.startcodon,
            self.stopcodon
        )

    def __repr__(self):
        return "%s(%d, %d, %r, %r, %r)" % (
            self.__class__.__name__,
            self.start,
            self.end,
            self.intact,
            self.startcodon,
            self.stopcodon
        )


class CDSCollection(object):
    def __init__(self, transcript, cdslist, part):
        self.transcript = transcript
        self.cdslist = cdslist
        self.part = part

    def __str__(self):
        return "%s(transcript_id = %s, transcript_name = %s, part= %s)" % (
            self.__class__.__name__,
            self.transcript.transcript_id,
            self.transcript.transcript_name,
            self.part
        )

    def __repr__(self):
        return "%s(%s, %s, %s)" % (
            self.__class__.__name__,
            self.transcript.transcript_id,
            self.transcript.transcript_name,
            self.part
        )

    @property
    def cds_range(self):
        if self.part == '5':
            start = 0
            end = start + self.cut_length - 1
        else:
            end = self.comp_length - 1
            # pyensembl coding_sequence_ranges does't not include the stop codon
            # but coding_sequence includes the stop codon, so we should minus 3
            # additionally for 3' part transcript
            start = end - self.cut_length + 1 - 3
        return start, end

    @property
    def transcript_id(self):
        """
        :return: transcript id of this collection of exons
        """
        return self.transcript.transcript_id

    @property
    def strand(self):
        return self.transcript.strand

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
    def nt_sequence(self):
        seq_start, seq_end = self.cds_range
        return self.transcript.coding_sequence[seq_start: seq_end+1]

    @property
    def cut_length(self):
        return sum([(cds.end-cds.start+1) for cds in self.cdslist])

    @property
    def comp_length(self):
        return len(self.transcript.coding_sequence)


class SVFusion(object):
    def __init__(self, sv, cdscollection_1, cdscollection_2):
        self.sv = sv
        self.cc_1 = cdscollection_1
        self.cc_2 = cdscollection_2
        self.nt_sequence = None
        self.aa_sequence = None
        self.neoepitopes = None

    def __str__(self):
        return "%s(sv = %r, cc_1 = %r, cc_2 = %r, nt_sequence = %s, aa_sequence = %s)" % (
            self.__class__.__name__,
            self.sv,
            self.cc_1,
            self.cc_2,
            self.nt_sequence,
            self.aa_sequence
        )

    def __repr__(self):
        return "%s(%r, %r, %r, %s, %s)" % (
            self.__class__.__name__,
            self.sv,
            self.cc_1,
            self.cc_2,
            self.nt_sequence,
            self.aa_sequence
        )

    def __eq__(self, other):
        return self.sv.sorted_coord == other.sv.sorted_coord

    def __hash__(self):
        return hash(self.sv.sorted_coord)

    def is_empty(self):
        if not self.cc_1 or not self.cc_2:
            return True

    @property
    def nt_sequence_ins(self):
        # the insertion sequence is for the forward strand,
        # so should be adjusted by the direction of final transcript
        if self.cc_1.strand == '+':
            return self.sv.insertion
        else:
            return str(Seq(self.sv.insertion).reverse_complement())

    @property
    def nt_sequence_cds(self):
        if self.cc_1.part == '5':
            return self.cc_1.nt_sequence + self.nt_sequence_ins + self.cc_2.nt_sequence
        else:
            return self.cc_2.nt_sequence + self.nt_sequence_ins + self.cc_1.nt_sequence

    @property
    def nt_sequence_3utr(self):
        if self.cc_1.part == '5':
            return self.cc_2.transcript.three_prime_utr_sequence
        else:
            return self.cc_1.transcript.three_prime_utr_sequence

    @property
    def frame_effect(self):
        # fixed a bug: the 5' cds collection could be empty, or less than 3 aa,
        # in this case, traditional start codon is lost and the program will search for the next start codon automatically.
        # but such prediction is of low reliability, we should annotate it and remove these fusions when necessary.
        if self.cc_1.part == '5' and not self.cc_1.nt_sequence.startswith('ATG'):
            return 'Start-loss'
        if self.cc_2.part == '5' and not self.cc_2.nt_sequence.startswith('ATG'):
            return 'Start-loss'
        if len(self.nt_sequence_cds) == 3*(len(self.aa_sequence)+1):
            return 'In-frame'
        elif len(self.nt_sequence_cds) > 3*(len(self.aa_sequence)+1):
            return 'Stop-gain'
        else:
            return 'Stop-loss'

    def output(self):
        tran_id_1 = self.cc_1.transcript_id
        gene_name_1 = self.cc_1.gene_name
        tran_id_2 = self.cc_2.transcript_id
        gene_name_2 = self.cc_2.gene_name
        chrom_1, chrom_2 = self.sv.chrom1, self.sv.chrom2
        pos_1, pos_2 = self.sv.pos1, self.sv.pos2
        return [chrom_1, str(pos_1), gene_name_1, tran_id_1,
                chrom_2, str(pos_2), gene_name_2, tran_id_2,
                str(self.sv.pattern), self.sv.svtype, self.frame_effect]