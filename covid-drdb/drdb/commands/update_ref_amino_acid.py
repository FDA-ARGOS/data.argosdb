import click
import requests
from pathlib import Path
from typing import Dict, Tuple, List, Any
from more_itertools import chunked

from ..cli import cli
from ..utils.csvv import dump_csv
from ..utils.codonutils import translate_codon

GENES_URL: str = (
    'https://raw.githubusercontent.com/hivdb/sierra-sars2/'
    'master/src/main/resources/genes.json'
)

NUCREF_URL: str = (
    'https://raw.githubusercontent.com/hivdb/sierra-sars2/'
    'master/src/main/resources/alignment-config.json'
)


def read_codons() -> Dict[Tuple[str, int], str]:
    refseq: str
    geneseq: str
    genedef: Dict[str, Any]
    fragment: Dict[str, Any]
    ref_start: int
    ref_end: int
    gene: str
    pos0: int
    bps: List[str]
    resp: requests.Response = requests.get(NUCREF_URL)
    config: Dict[str, Any] = resp.json()
    genedefs: List[Dict[str, Any]] = []
    result: Dict[Tuple[str, int], str] = {}
    for fragment in config['fragmentConfig']:
        if fragment['fragmentName'] == 'Wuhan-Hu-1::NC_045512.2':
            refseq = fragment['refSequence']
        elif 'geneName' in fragment:
            genedefs.append(fragment)
    for genedef in genedefs:
        gene = genedef['geneName'][5:]
        geneseq = ''.join([
            refseq[ref_start - 1:ref_end]
            for ref_start, ref_end in genedef['refRanges']
        ])
        for pos0, bps in enumerate(chunked(geneseq, 3)):
            result[(gene, pos0 + 1)] = ''.join(bps)
    return result


@cli.command()
@click.argument(
    'payload_dir',
    type=click.Path(
        dir_okay=True,
        exists=True,
        file_okay=False
    )
)
def update_ref_amino_acid(payload_dir: str) -> None:
    refaa_csv: Path = Path(payload_dir) / 'tables' / 'ref_amino_acid.csv'
    resp: requests.Response = requests.get(GENES_URL)
    genes: List[Dict[str, Any]] = resp.json()
    rows: List[Dict[str, Any]] = []
    codons: Dict[Tuple[str, int], str] = read_codons()
    for gene_def in genes:
        gene: str = gene_def['abstractGene']
        refseq: str = gene_def['refSequence']
        for pos0, aa in enumerate(refseq):
            pos = pos0 + 1
            codon = codons[(gene, pos)]
            aa_from_codon = translate_codon(
                codon.encode('ASCII')
            ).decode('ASCII')
            if aa != aa_from_codon:
                click.echo(
                    f"Error: {gene} position {pos}'s ref amino acid "
                    f"{aa} does not match its codon {codon} ({aa_from_codon})",
                    err=True
                )
                raise click.Abort()
            rows.append({
                'gene': gene,
                'position': pos,
                'amino_acid': aa,
                'codon': codons[(gene, pos0 + 1)]
            })
    dump_csv(
        refaa_csv,
        rows,
        ['gene', 'position', 'amino_acid', 'codon']
    )
