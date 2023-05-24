import click
from pathlib import Path
from ..cli import cli
from ..utils.csvv import load_csv, dump_csv


def merge_expgroup_core(tables_dir):
    assays = tables_dir / 'assays.csv'
    rxpots = tables_dir / 'rx_potency'
    expgrps = tables_dir / 'experiment_groups'
    vtype_lookup = {
        a['assay_name']: a['virus_type']
        for a in load_csv(assays)
    }
    lookup = {}
    for expgrp in expgrps.iterdir():
        if expgrp.suffix.lower() != '.csv':
            click.echo('Skip {}'.format(expgrp))
            continue
        rows = load_csv(expgrp)
        for row in rows:
            lookup[(
                row['ref_name'],
                row['virus_type'],
                row['potency_type']
            )] = row
    for rxpot in rxpots.iterdir():
        if rxpot.suffix.lower() != '.csv':
            click.echo('Skip {}'.format(rxpot))
            continue
        rows = load_csv(rxpot)
        for row in rows:
            key = (
                row['ref_name'],
                vtype_lookup[row['assay_name']],
                row['potency_type']
            )
            if key in lookup:
                expgrp_row = lookup[key]
                row['potency_upper_limit'] = expgrp_row['potency_upper_limit']
                row['potency_lower_limit'] = expgrp_row['potency_lower_limit']
                row['potency_unit'] = expgrp_row['potency_unit']
        dump_csv(
            rxpot,
            records=rows,
            headers=[
                'ref_name',
                'rx_name',
                'iso_name',
                'section',
                'assay_name',
                'potency_type',
                'potency',
                'cumulative_count',
                'potency_upper_limit',
                'potency_lower_limit',
                'potency_unit',
                'date_added',
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
def merge_expgroup(payload_dir):
    payload_dir = Path(payload_dir)
    merge_expgroup_core(payload_dir / 'tables')
