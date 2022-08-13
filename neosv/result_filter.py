import os
class StructuralVariant(object):
    """
    Class for storing SV information:
    1. chromosomes, format: string
    2. positions, format: integer
    2. insertion sequence, format: string
    3. sv pattern, format: integer
    NOTE: sv pattern corresponds to 4 patterns in VCF specification
    """
    def __init__(self, chrom1, pos1, chrom2, pos2, info):
        self.chrom1 = str(chrom1).replace('chr', '')
        self.pos1 = int(pos1)
        self.chrom2 = str(chrom2).replace('chr', '')
        self.pos2 = int(pos2)
        self.info = info

    def __str__(self):
        return "%s(chrom1 = %s, pos1 = %d, chrom2 = %s, pos2 = %d)" % (
            self.__class__.__name__,
            self.chrom1,
            self.pos1,
            self.chrom2,
            self.pos2
        )

    def __repr__(self):
        return "%s(%s, %d, %s, %d)" % (
            self.__class__.__name__,
            self.chrom1,
            self.pos1,
            self.chrom2,
            self.pos2
        )

    def __eq__(self, other):
        if self.bp1_sorted[0] != other.bp1_sorted[0] or self.bp2_sorted[0] != other.bp2_sorted[0]:
            return False
        else:
            if (other.bp1_sorted[1]-100) <= self.bp1_sorted[1] <= (other.bp1_sorted[1]+100) and (other.bp2_sorted[1]-100) <= self.bp2_sorted[1] <= (other.bp2_sorted[1]+100):
                return True
            else:
                return False

    @property
    def sorted_coord(self):
        """
        :return: sorted coordinates of two breakpoints, used for class identity
        """
        bp1 = [self.chrom1, self.pos1]
        bp2 = [self.chrom2, self.pos2]
        return sorted([bp1, bp2], key = lambda bp: (bp[0], bp[1]))

    @property
    def bp1_sorted(self):
        return self.sorted_coord[0]

    @property
    def bp2_sorted(self):
        return self.sorted_coord[1]


def anno_filter(annofile):
    svlist = []
    with open(annofile, 'r') as f:
        next(f)
        for line in f:
            tmpline = line.rstrip().split('\t')
            chrom1 = tmpline[0]
            pos1 = tmpline[1]
            chrom2 = tmpline[7]
            pos2 = tmpline[8]
            sv = StructuralVariant(chrom1, pos1, chrom2, pos2, line)
            if sv not in svlist:
                svlist.append(sv)
    return svlist


def neo_filter(neofile, annolist):
    sv_coord_list = [sv.chrom1 + '_' + str(sv.pos1) + '_' + sv.chrom2 + '_' + str(sv.pos2) for sv in annolist]
    svlist = []
    with open(neofile, 'r') as f:
        next(f)
        for line in f:
            tmpline = line.rstrip().split('\t')
            chrom1 = tmpline[0]
            pos1 = tmpline[1]
            chrom2 = tmpline[4]
            pos2 = tmpline[5]
            sv_coord = chrom1 + '_' + pos1 + '_' + chrom2 + '_' + pos2
            sv = StructuralVariant(chrom1, pos1, chrom2, pos2, line)
            if sv_coord in sv_coord_list:
                svlist.append(sv)
    return svlist


def write_neofile(svlist, filepath):
    header = '\t'.join(['chrom1', 'pos1', 'gene1', 'transcript_id1',
                            'chrom2', 'pos2', 'gene2', 'transcript_id2',
                            'svpattern', 'svtype', 'frameshift',
                            'neoantigen', 'allele', 'affinity', 'rank'])
    with open(filepath, 'w') as f:
        f.write(header + '\n')
        for sv in svlist:
            f.write(sv.info)


def write_annofile(svlist, filepath):
    header = '\t'.join(['chrom1', 'pos1', 'function1', 'gene1', 'transcript_id1', 'strand1', 'transcript_retain1',
                        'chrom2', 'pos2', 'function2', 'gene2', 'transcript_id2', 'strand2', 'transcript_retain2',
                        'svpattern', 'svtype', 'fusion'])
    with open(filepath, 'w') as f:
        f.write(header + '\n')
        for sv in svlist:
            f.write(sv.info)


indir = '/home/ruibinxi_pkuhpc/lustre1/shiyang/SV_neoantigen/svanno_result/'

samples = [f for f in os.listdir(indir)]
for sample in samples:
    print(sample)
    file_anno = os.path.join(indir, sample, sample + '.anno.txt')
    file_neo = os.path.join(indir, sample, sample + '.neoantigen.txt')
    annolist = anno_filter(file_anno)
    neolist = neo_filter(file_neo, annolist)
    file_anno_filtered = os.path.join(indir, sample, sample + '.anno.filtered.txt')
    file_neo_filtered = os.path.join(indir, sample, sample + '.neoantigen.filtered.txt')
    write_annofile(annolist, file_anno_filtered)
    write_neofile(neolist, file_neo_filtered)
