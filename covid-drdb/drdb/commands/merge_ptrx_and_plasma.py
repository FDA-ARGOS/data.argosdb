import click
from pathlib import Path
from ..cli import cli
from ..utils.csvv import load_csv, dump_csv


def merge_ptrx_and_plasma_core(tables_dir):
    ptrx_dir = tables_dir / 'patient_treatments'
    rx_plasma_dir = tables_dir / 'rx_plasma'
    lookup = {}
    for ptrx in ptrx_dir.iterdir():
        if ptrx.suffix.lower() != '.csv':
            click.echo('Skip {}'.format(ptrx))
            continue
        rows = load_csv(ptrx)
        for row in rows:
            lookup[(row['ref_name'], row['rx_name'])] = row['patient_name']
    for rxp in rx_plasma_dir.iterdir():
        if rxp.suffix.lower() != '.csv':
            click.echo('Skip {}'.format(rxp))
            continue
        rows = load_csv(rxp)
        for row in rows:
            key = (row['ref_name'], row['rx_name'])
            row['patient_name'] = lookup.get(key)
        dump_csv(
            rxp,
            records=rows,
            headers=[
                'ref_name',
                'rx_name',
                'patient_name',
                'titer',
                'collection_date',
                'cumulative_group'
            ],
            BOM=True
        )


@cli.command()
@click.argument(
    'payload_dir',
    type=click.Path(
        dir_okay=True,
        exists=True,
        file_okay=False
    )
)
def merge_ptrx_and_plasma(payload_dir):
    payload_dir = Path(payload_dir)
    merge_ptrx_and_plasma_core(payload_dir / 'tables')
    merge_ptrx_and_plasma_core(payload_dir / 'excluded')
