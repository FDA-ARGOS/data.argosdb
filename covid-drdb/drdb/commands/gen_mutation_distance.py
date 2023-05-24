import click
from pathlib import Path
from itertools import product
from typing import List, Iterable, Dict, Any

from ..cli import cli
from ..utils.csvv import dump_csv
from ..utils.codonutils import (
    AMINO_ACID_LIST,
    REVERSE_CODON_TABLE
)


def get_bp_changes(
    source_codon: Iterable[int],
    dest_codon: Iterable[int]
) -> int:
    a: int
    b: int
    ndiffs: int = 0
    for a, b in zip(source_codon, dest_codon):
        ndiffs += a != b
    return ndiffs


@cli.command()
@click.argument(
    'payload_dir',
    type=click.Path(
        dir_okay=True,
        exists=True,
        file_okay=False
    )
)
def gen_mutation_distance(payload_dir: str) -> None:
    codon: bytes
    dest_codon: bytes
    dest_codons: List[bytes]
    source_aa: int
    aa: int
    rows: List[Dict[str, Any]] = []
    mutdist_csv: Path = Path(payload_dir) / 'tables' / 'mutation_distance.csv'
    codons: List[bytes] = [bytes(cd) for cd in product(b'ACGT', repeat=3)]
    for codon, aa in product(codons, AMINO_ACID_LIST):
        dest_codons = REVERSE_CODON_TABLE[aa]
        bp_changes = sorted([
            get_bp_changes(codon, dest_codon)
            for dest_codon in dest_codons
        ])
        rows.append({
            'source_codon': codon.decode('ASCII'),
            'dest_amino_acid': 'stop' if aa == 42 else chr(aa),
            'min_bp_changes': bp_changes[0],
            'max_bp_changes': bp_changes[-1]
        })
    dump_csv(
        mutdist_csv,
        rows,
        ['source_codon',
         'dest_amino_acid',
         'min_bp_changes',
         'max_bp_changes']
    )
