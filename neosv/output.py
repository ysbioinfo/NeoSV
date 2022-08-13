def write_annot(filepath, sveffects):
    """
    :param filepath: outfile for annotation
    :param sveffects: a list of sveffect classes
    :return: None
    """
    header = '\t'.join(['chrom1', 'pos1', 'function1', 'gene1', 'transcript_id1', 'strand1', 'transcript_retain1',
                        'chrom2', 'pos2', 'function2', 'gene2', 'transcript_id2', 'strand2', 'transcript_retain2',
                        'svpattern', 'svtype', 'fusion'])
    with open(filepath, 'w') as f:
        f.write(header + '\n')
        for sveffect in sveffects:
            f.write('\t'.join(sveffect.output()) + '\n')


def write_fusion(filepath, svfusions, dict_neo):
    """
    :param filepath: outfile for fusion derived neoantigen
    :param svfusions: a list of svfusion classes
    :param dict_neo: a dictionary for neoepitopes
    :return: None
    """
    with open(filepath, 'w') as f:
        header = '\t'.join(['chrom1', 'pos1', 'gene1', 'transcript_id1',
                            'chrom2', 'pos2', 'gene2', 'transcript_id2',
                            'svpattern', 'svtype', 'frameshift',
                            'neoantigen', 'allele', 'affinity', 'rank'])
        f.write(header + '\n')
        for svfusion in svfusions:
            for neoepitope in svfusion.neoepitopes:
                for allele in dict_neo[neoepitope]:
                    if dict_neo[neoepitope][allele][2] == 'PASS':
                        affinity = str(dict_neo[neoepitope][allele][0])
                        rank = str(dict_neo[neoepitope][allele][1])
                        f.write('\t'.join(['\t'.join(svfusion.output()), neoepitope, allele, affinity, rank]) + '\n')