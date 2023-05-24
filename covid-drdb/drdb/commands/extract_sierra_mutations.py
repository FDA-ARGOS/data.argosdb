import os
import csv
import click
from tqdm import tqdm  # type: ignore
from typing import TextIO, Dict, List, Iterable
from itertools import groupby

from ..cli import cli

HEADER: List[str] = [
    'iso_name',
    'gene',
    'position',
    'amino_acid',
    'count',
    'total'
]


def remove_mixtures(rows: Iterable[Dict[str, str]]) -> List[Dict[str, str]]:
    results: List[Dict[str, str]] = []
    rows = list(rows)
    has_percent = any('Percent' in r for r in rows)
    for (gene, pos), part in groupby(
        rows, lambda r: (r['Gene'], r['Position'])
    ):
        partlst = [
            r for r in part
            if r['RefAA'] != r['MutAA'] and
            r['MutAA'] != 'X']
        if has_percent:
            partlst = [r for r in partlst if float(r['Percent']) > 20]
        elif len(partlst) > 1:
            continue
        results.extend(partlst)
    return results


@cli.command()
@click.argument(
    'sierra-mutation-list-dir',
    type=click.Path(exists=True, file_okay=False)
)
@click.argument(
    'output_mutations_csv',
    type=click.File('w', encoding='UTF-8-sig')
)
def extract_sierra_mutations(
    sierra_mutation_list_dir: str,
    output_mutations_csv: TextIO
) -> None:
    fp: TextIO
    mutlist: str
    rows: List[Dict[str, str]]
    row: Dict[str, str]

    reader: csv.DictReader
    writer: csv.DictWriter = csv.DictWriter(
        output_mutations_csv,
        HEADER
    )
    writer.writeheader()

    for mutlist in tqdm(os.listdir(sierra_mutation_list_dir)):
        if not mutlist.startswith('MutationList') or \
                not mutlist.endswith('.csv'):
            continue
        mutlist = os.path.join(sierra_mutation_list_dir, mutlist)
        with open(mutlist, encoding='UTF-8-sig') as fp:
            reader = csv.DictReader(fp)
            rows = remove_mixtures(reader)
            for row in rows:
                writer.writerow({
                    'iso_name': row['Sequence Name'],
                    'gene': row['Gene'],
                    'position': row['Position'],
                    'amino_acid': row['MutAA'],
                    'count': 'NULL',
                    'total': 'NULL'
                })
