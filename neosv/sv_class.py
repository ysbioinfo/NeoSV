class StructuralVariant(object):
    """
    Class for storing SV information:
    1. chromosomes, format: string
    2. positions, format: integer
    2. insertion sequence, format: string
    3. sv pattern, format: integer
    NOTE: sv pattern corresponds to 4 patterns in VCF specification
    """
    def __init__(self, chrom1, pos1, chrom2, pos2, insertion, pattern):
        self.chrom1 = str(chrom1).replace('chr', '')
        self.pos1 = int(pos1)
        self.chrom2 = str(chrom2).replace('chr', '')
        self.pos2 = int(pos2)
        self.insertion = insertion
        self.pattern = int(pattern)

    def __str__(self):
        return "%s(chrom1 = %s, pos1 = %d, chrom2 = %s, pos2 = %d, insertion = '%s', pattern = %d)" % (
            self.__class__.__name__,
            self.chrom1,
            self.pos1,
            self.chrom2,
            self.pos2,
            self.insertion,
            self.pattern
        )

    def __repr__(self):
        return "%s(%s, %d, %s, %d, '%s', %d)" % (
            self.__class__.__name__,
            self.chrom1,
            self.pos1,
            self.chrom2,
            self.pos2,
            self.insertion,
            self.pattern
        )

    def __eq__(self, other):
        return self.sorted_coord == other.sorted_coord

    def __hash__(self):
        return hash(self.sorted_coord)

    @property
    def sorted_coord(self):
        """
        :return: sorted coordinates of two breakpoints, used for class identity
        """
        bp1 = self.chrom1 + '_' + str(self.pos1)
        bp2 = self.chrom2 + '_' + str(self.pos2)
        return tuple(sorted([bp1, bp2]))

    @property
    def svtype(self):
        if self.chrom1 != self.chrom2:
            return 'TRA'
        else:
            if self.pos1 < self.pos2:
                if self.pattern == 1:
                    return 'DEL'
                elif self.pattern == 2:
                    return 'h2hINV'
                elif self.pattern == 3:
                    return 'DUP'
                else:
                    return 't2tINV'
            else:
                if self.pattern == 1:
                    return 'DUP'
                elif self.pattern == 2:
                    return 'h2hINV'
                elif self.pattern == 3:
                    return 'DEL'
                else:
                    return 't2tINV'
