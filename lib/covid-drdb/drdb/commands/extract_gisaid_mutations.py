import re
import csv
import click
from tqdm import tqdm  # type: ignore
from typing import TextIO, Dict, List, NamedTuple, Optional, Set

from ..cli import cli

GISAID_GENES: Dict[str, str] = {
    'NSP1': 'nsp1',
    'NSP2': 'nsp2',
    'NSP3': 'PLpro',
    'NSP4': 'nsp4',
    'NSP5': '_3CLpro',
    'NSP6': 'nsp6',
    'NSP7': 'nsp7',
    'NSP8': 'nsp8',
    'NSP9': 'nsp9',
    'NSP10': 'nsp10',
    'NSP11': 'RdRP',
    'NSP12': 'RdRP',
    'NSP13': 'nsp13',
    'NSP14': 'nsp14',
    'NSP15': 'nsp15',
    'NSP16': 'nsp16',
    'Spike': 'S',
    'NS3': 'ORF3a',
    'E': 'E',
    'M': 'M',
    'NS6': 'ORF6',
    'NS7a': 'ORF7a',
    'NS7b': 'ORF7b',
    'NS8': 'ORF8',
    'N': 'N',
    'NS10': 'ORF10'
}

ORDERED_GENES: List[str] = [
    'nsp1', 'nsp2', 'PLpro', 'nsp4', '_3CLpro', 'nsp5', 'nsp6',
    'nsp7', 'nsp8', 'nsp9', 'nsp10', 'RdRP', 'nsp13', 'nsp14',
    'nsp15', 'nsp16', 'S', 'ORF3a', 'E', 'M', 'ORF6', 'ORF7a',
    'ORF7b', 'ORF8', 'N', 'ORF10'
]

HEADER: List[str] = [
    'Accession',
    'Gene',
    'Position',
    'MutAA'
]


class Mutation(NamedTuple):
    gene: str
    position: int
    amino_acid: str


def parse_mutations(mutations: str) -> List[Mutation]:
    gene: str
    pos: str
    refaa: str
    mutaa: str
    muttext: str
    match: Optional[re.Match]

    mutlist: List[Mutation] = []
    mutations = mutations.strip('()')
    for muttext in mutations.split(','):
        if not muttext:
            continue
        match = re.match(r'([A-Za-z\d]+)_([A-Za-z]+)(\d+)([A-Za-z]+)', muttext)
        if not match:
            click.echo(f'Invalid mutation: {muttext}', err=True)
            continue
        gene, refaa, pos, mutaa = match.groups()
        if refaa == 'ins':
            mutaa = 'ins'
        mutlist.append(Mutation(GISAID_GENES[gene], int(pos), mutaa))
    return sorted(mutlist, key=lambda mut: (
        ORDERED_GENES.index(mut.gene),
        mut.position,
        mut.amino_acid
    ))


@cli.command()
@click.argument(
    'input_metadata_tsv',
    type=click.File(encoding='UTF-8-sig')
)
@click.argument(
    'output_mutations_csv',
    type=click.File('w', encoding='UTF-8-sig')
)
@click.option(
    '--gene',
    type=str,
    help=(
        'Filter gene(s) if specified in '
        'no-space, comma-separated format.'
    )
)
def extract_gisaid_mutations(
    input_metadata_tsv: TextIO,
    output_mutations_csv: TextIO,
    gene: str
) -> None:
    row: Dict[str, str]
    mut: Mutation
    mutations: List[Mutation]
    genes: Set = set(gene.split(','))
    for gene in genes:
        if gene not in ORDERED_GENES:
            raise click.BadOptionUsage(
                option_name='gene',
                message=(
                    'Invalid gene specified: {}. Valid genes: {}'
                    .format(gene, ','.join(ORDERED_GENES))
                )
            )

    reader: csv.DictReader = csv.DictReader(input_metadata_tsv, delimiter='\t')
    writer: csv.DictWriter = csv.DictWriter(
        output_mutations_csv,
        HEADER
    )
    writer.writeheader()
    for row in tqdm(reader, unit='seqs'):
        mutations = parse_mutations(row['AA Substitutions'])
        for mut in mutations:
            if genes and mut.gene not in genes:
                continue
            writer.writerow({
                'Accession': row['Accession ID'],
                'Gene': mut.gene,
                'Position': mut.position,
                'MutAA': mut.amino_acid
            })
