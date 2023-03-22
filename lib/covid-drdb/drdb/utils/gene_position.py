import warnings
from typing import Tuple


def translate_gene_position(gene: str, pos: int) -> Tuple[str, int]:
    lower_gene = gene.lower()
    if lower_gene in ('orf1a', 'orf1ab'):
        if pos <= 180:
            return 'nsp1', pos
        elif pos <= 818:
            return 'nsp2', pos - 180
        elif pos <= 2763:
            return 'PLpro', pos - 818
        elif pos <= 3263:
            return 'nsp4', pos - 2763
        elif pos <= 3569:
            return '_3CLpro', pos - 3263
        elif pos <= 3859:
            return 'nsp6', pos - 3569
        elif pos <= 3942:
            return 'nsp7', pos - 3859
        elif pos <= 4140:
            return 'nsp8', pos - 3942
        elif pos <= 4253:
            return 'nsp9', pos - 4140
        elif pos <= 4392:
            return 'nsp10', pos - 4253
        else:
            return 'RdRP', pos - 4392
    elif lower_gene == 'nsp3':
        return 'PLpro', pos
    elif lower_gene in ('nsp5', '3clpro', '3cl', 'mpro', 'mainpro'):
        return '_3CLpro', pos
    elif lower_gene == 'nsp11':
        return 'RdRP', pos
    elif lower_gene == 'nsp12':
        warnings.warn(
            'Unchanged ambiguous gene nsp12: '
            'nsp12 can be either the whole RdRP '
            'or just partial RdRP in ORF1b.'
        )
        return 'nsp12', pos
    elif lower_gene.startswith('nsp'):
        return lower_gene, pos
    elif lower_gene == 'orf1b':
        if pos <= 923:
            return 'RdRP', pos + 9
        elif pos <= 1524:
            return 'nsp13', pos - 923
        elif pos <= 2051:
            return 'nsp14', pos - 1524
        elif pos <= 2397:
            return 'nsp15', pos - 2051
        else:
            return 'nsp16', pos - 2397
    elif lower_gene in ('s', 'spike', 'ns2', 'orf2', 'gp02'):
        return 'S', pos
    elif lower_gene in ('ns3', 'orf3', 'gp03'):
        return 'ORF3', pos
    elif lower_gene in ('e', 'ns4', 'orf4', 'gp04'):
        return 'E', pos
    elif lower_gene in ('m', 'ns5', 'orf5', 'gp05'):
        return 'M', pos
    elif lower_gene in ('ns6', 'orf6', 'gp06'):
        return 'ORF6', pos
    elif lower_gene in ('ns7a', 'orf7a', 'gp07'):
        return 'ORF7a', pos
    elif lower_gene in ('ns7b', 'orf7b', 'gp08'):
        return 'ORF7b', pos
    elif lower_gene in ('ns8', 'orf8', 'gp09'):
        return 'ORF8', pos
    elif lower_gene in ('n', 'ns9', 'orf9', 'gp10'):
        return 'N', pos
    elif lower_gene in ('ns10', 'orf10', 'gp11'):
        return 'ORF10', pos
    return gene, pos
