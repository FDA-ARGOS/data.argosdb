import click
from pathlib import Path
from ..cli import cli
from ..utils.csvv import load_csv, dump_csv


def merge_cpvp_core(tables_dir):
    rxcps = tables_dir / 'rx_conv_plasma'
    rxvps = tables_dir / 'rx_vacc_plasma'
    rx_plasma_dir = tables_dir / 'rx_plasma'
    merged = {}
    for rxcp in rxcps.iterdir():
        if rxcp.suffix.lower() != '.csv':
            click.echo('Skip {}'.format(rxcp))
            continue
        rows = load_csv(rxcp)
        merged[rows[0]['ref_name']] = rows
    for rxvp in rxvps.iterdir():
        if rxvp.suffix.lower() != '.csv':
            click.echo('Skip {}'.format(rxvp))
            continue
        rows = load_csv(rxvp)
        ref_name = rows[0]['ref_name']
        if ref_name in merged:
            merged[ref_name] = merged[ref_name] + rows
        else:
            merged[ref_name] = rows
    for ref_name, rows in merged.items():
        rx_plasma = rx_plasma_dir / '{}-plasma.csv'.format(ref_name.lower())
        dump_csv(
            rx_plasma,
            records=rows,
            headers=[
                'ref_name',
                'rx_name',
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
def merge_cpvp(payload_dir):
    payload_dir = Path(payload_dir)
    merge_cpvp_core(payload_dir / 'tables')
    merge_cpvp_core(payload_dir / 'excluded')
