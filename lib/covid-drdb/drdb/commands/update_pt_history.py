import click
from pathlib import Path
from itertools import chain
from ..cli import cli
from ..utils.csvv import load_csv, dump_csv, load_multiple_csvs


def load_rx2pt(tables_dir):
    ptrx = load_multiple_csvs(tables_dir / 'patient_treatments')
    return {
        (row['ref_name'], row['rx_name']): row['patient_name']
        for row in ptrx
    }


def update_pt_history_core(payload_dir):
    tables_dir = payload_dir / 'tables'
    excluded_dir = payload_dir / 'excluded'
    rxvp_dir = tables_dir / 'rx_vacc_plasma'
    rxvp_dir_exc = excluded_dir / 'rx_vacc_plasma'
    pthistory_dir = tables_dir / 'patient_history'
    rx2pt = load_rx2pt(tables_dir)
    vaccine_name_lookup = {}
    found = set()
    for rxvp in chain(rxvp_dir.iterdir(), rxvp_dir_exc.iterdir()):
        if rxvp.suffix.lower() != '.csv':
            click.echo('Skip {}'.format(rxvp))
            continue
        for row in load_csv(rxvp):
            rxkey = (row['ref_name'], row['rx_name'])
            if row['vaccine_name']:
                if rxkey not in rx2pt:
                    click.echo('No patient exists: {}'.format(rxkey))
                    continue
                key = (
                    row['ref_name'],
                    rx2pt[(row['ref_name'], row['rx_name'])]
                )
                vacc = row['vaccine_name']
                if (
                    key in vaccine_name_lookup and
                    vaccine_name_lookup[key] != vacc
                ):
                    raise ValueError(
                        'Different vaccine name for same patient: {} ({})'
                        .format(key, vacc)
                    )
                vaccine_name_lookup[key] = vacc
            else:
                click.echo('No vaccine name in rxVP: {}'.format(rxkey))
    for pth in pthistory_dir.iterdir():
        if pth.suffix.lower() != '.csv':
            click.echo('Skip {}'.format(pth))
            continue
        rows = load_csv(pth)
        for row in rows:
            if row['event'] in (
                "1st dose",
                "1st dose isolation",
                "2nd dose",
                "2nd dose isolation",
                "2nd dose saliva isolation",
                "3rd dose",
                "3rd dose isolation"
            ):
                key = (row['ref_name'], row['patient_name'])
                row['vaccine_name'] = vaccine_name_lookup[key]
                found.add(key)
        dump_csv(
            pth,
            records=rows,
            headers=[
                'ref_name',
                'patient_name',
                'event',
                'event_date_cmp',
                'event_date',
                'location',
                'iso_name',
                'vaccine_name',
                'severity'
            ]
        )
    click.echo('Missing patient: {}'.format(
        set(vaccine_name_lookup.keys()) - found
    ))


@cli.command()
@click.argument(
    'payload_dir',
    type=click.Path(
        dir_okay=True,
        exists=True,
        file_okay=False
    )
)
def update_pt_history(payload_dir):
    payload_dir = Path(payload_dir)
    update_pt_history_core(payload_dir)
