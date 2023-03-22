import click
from pathlib import Path
from ..cli import cli
from ..utils.csvv import load_csv, dump_csv, load_multiple_csvs

from itertools import chain


def load_isolate_variant_lookup(tables_dir):
    isolates = load_csv(tables_dir / 'isolates.csv')
    isolates += load_multiple_csvs(tables_dir / 'isolates.d')
    return {
        iso['iso_name']: iso['var_name']
        for iso in isolates
        if iso['var_name']
    }


def load_rx_plasma_lookup(tables_dir):
    rx_plasma = load_multiple_csvs(tables_dir / 'rx_plasma')
    lookup = {}
    for rxp in rx_plasma:
        key = (
            rxp['ref_name'],
            rxp['subject_name'],
            rxp['collection_date'])
        if key not in lookup:
            lookup[key] = []
        lookup[key].append(rxp)
    return lookup


def pth_to_infection(row, var_lookup):
    return {
        'ref_name': row['ref_name'],
        'subject_name': row['subject_name'],
        'infection_date_cmp': row['event_date_cmp'],
        'infection_date': row['event_date'],
        'infected_var_name': (
            var_lookup.get(row['iso_name']) or row['iso_name']
        ),
        'location': row['location'],
        'section': None
    }


def pth_to_severity(row):
    return {
        'ref_name': row['ref_name'],
        'subject_name': row['subject_name'],
        'start_date_cmp': row['event_date_cmp'],
        'start_date': row['event_date'],
        'end_date_cmp': None,
        'end_date': None,
        'severity': row['severity'],
        'section': None
    }


def pth_to_plasma(row, rx_plasma_lookup):
    key = (
        row['ref_name'],
        row['subject_name'],
        row['event_date']
    )
    rxps = rx_plasma_lookup.pop(key, None)
    if not rxps:
        click.echo('rx_plasma not exist for {}'.format(key), err=True)
        yield {
            'ref_name': row['ref_name'],
            'subject_name': row['subject_name'],
            'rx_name': None,
            'collection_date_cmp': row['event_date_cmp'],
            'collection_date': row['event_date'],
            'location': row['location'],
            'cumulative_group': None,
            'section': None
        }
    else:
        for rxp in rxps:
            yield {
                'ref_name': row['ref_name'],
                'subject_name': row['subject_name'],
                'rx_name': rxp['rx_name'],
                'collection_date_cmp': row['event_date_cmp'],
                'collection_date': row['event_date'],
                'location': row['location'],
                'cumulative_group': rxp['cumulative_group'],
                'section': None
            }


def pth_to_isolate(row):
    return {
        'ref_name': row['ref_name'],
        'subject_name': row['subject_name'],
        'collection_date_cmp': row['event_date_cmp'],
        'collection_date': row['event_date'],
        'iso_name': row['iso_name'],
        'iso_source': None,
        'iso_culture': None,
        'location': row['location'],
        'section': None
    }


def pth_to_vaccine(row):
    return {
        'ref_name': row['ref_name'],
        'subject_name': row['subject_name'],
        'vaccination_date_cmp': row['event_date_cmp'],
        'vaccination_date': row['event_date'],
        'vaccine_name': row['vaccine_name'],
        'dosage': int(row['event'][0]),
        'location': row['location'],
        'section': None
    }


def refactor_sbj_history_core(tables_dir, var_lookup):
    pth_dir = tables_dir / 'subject_history'

    sbjinf_dir = tables_dir / 'subject_infections'
    sbjinf_dir.mkdir(exist_ok=True)

    sbjsev_dir = tables_dir / 'subject_severity'
    sbjsev_dir.mkdir(exist_ok=True)

    sbjp_dir = tables_dir / 'subject_plasma'
    sbjp_dir.mkdir(exist_ok=True)

    sbjiso_dir = tables_dir / 'subject_isolates'
    sbjiso_dir.mkdir(exist_ok=True)

    sbjvacc_dir = tables_dir / 'subject_vaccines'
    sbjvacc_dir.mkdir(exist_ok=True)

    rx_plasma_lookup = load_rx_plasma_lookup(tables_dir)

    for pth_file in pth_dir.iterdir():
        if pth_file.suffix.lower() != '.csv':
            click.echo('Skip {}'.format(pth_file))
            continue
        pth_rows = load_csv(pth_file)
        inf_rows = []
        sev_rows = {}
        plasma_rows = []
        iso_rows = []
        vacc_rows = []

        for row in pth_rows:
            sbjkey = (row['ref_name'], row['subject_name'])
            if row['event'] == 'infection':
                inf_rows.append(pth_to_infection(row, var_lookup))
                if row['severity']:
                    if sbjkey not in sev_rows:
                        sev_rows[sbjkey] = []
                    sev_rows[sbjkey].append(pth_to_severity(row))

            elif (
                (row['event'] == 'isolation' and row['iso_name'] is None) or
                row['event'].endswith('dose isolation')
            ):
                plasma_rows.extend(pth_to_plasma(row, rx_plasma_lookup))
                if sbjkey in sev_rows:
                    sev_rows[sbjkey][-1][
                        'end_date_cmp'] = row['event_date_cmp']
                    sev_rows[sbjkey][-1]['end_date'] = row['event_date']

            elif row['event'] == 'isolation' and row['iso_name']:
                iso_rows.append(pth_to_isolate(row))
                if sbjkey in sev_rows:
                    sev_rows[sbjkey][-1][
                        'end_date_cmp'] = row['event_date_cmp']
                    sev_rows[sbjkey][-1]['end_date'] = row['event_date']

            elif row['event'].endswith('dose'):
                vacc_rows.append(pth_to_vaccine(row))
                if sbjkey in sev_rows:
                    sev_rows[sbjkey][-1][
                        'end_date_cmp'] = row['event_date_cmp']
                    sev_rows[sbjkey][-1]['end_date'] = row['event_date']

        lower_ref_name = pth_rows[0]['ref_name'].lower()

        sbjinf_file = sbjinf_dir / '{}-inf.csv'.format(lower_ref_name)
        dump_csv(
            sbjinf_file,
            records=inf_rows,
            headers=[
                'ref_name',
                'subject_name',
                'infection_date_cmp',
                'infection_date',
                'infected_var_name',
                'location',
                'section'
            ],
            BOM=True
        )

        sbjsev_file = sbjsev_dir / '{}-sev.csv'.format(lower_ref_name)
        dump_csv(
            sbjsev_file,
            records=list(chain(*sev_rows.values())),
            headers=[
                'ref_name',
                'subject_name',
                'start_date_cmp',
                'start_date',
                'end_date_cmp',
                'end_date',
                'severity',
                'section'
            ],
            BOM=True
        )

        sbjp_file = sbjp_dir / '{}-plasma.csv'.format(lower_ref_name)
        dump_csv(
            sbjp_file,
            records=plasma_rows,
            headers=[
                'ref_name',
                'subject_name',
                'rx_name',
                'collection_date_cmp',
                'collection_date',
                'location',
                'cumulative_group',
                'section'
            ],
            BOM=True
        )

        sbjiso_file = sbjiso_dir / '{}-sbjiso.csv'.format(lower_ref_name)
        dump_csv(
            sbjiso_file,
            records=iso_rows,
            headers=[
                'ref_name',
                'subject_name',
                'collection_date_cmp',
                'collection_date',
                'iso_name',
                'iso_source',
                'iso_culture',
                'location',
                'section'
            ],
            BOM=True
        )

        sbjvacc_file = sbjvacc_dir / '{}-vacc.csv'.format(lower_ref_name)
        dump_csv(
            sbjvacc_file,
            records=vacc_rows,
            headers=[
                'ref_name',
                'subject_name',
                'vaccination_date_cmp',
                'vaccination_date',
                'vaccine_name',
                'dosage',
                'location',
                'section'
            ],
            BOM=True
        )
    if rx_plasma_lookup:
        click.echo('Unable to import following rx_plasma:', err=True)
        for key in rx_plasma_lookup:
            click.echo('  {}'.format(key), err=True)


@cli.command()
@click.argument(
    'payload_dir',
    type=click.Path(
        dir_okay=True,
        exists=True,
        file_okay=False
    )
)
def refactor_sbj_history(payload_dir):
    payload_dir = Path(payload_dir)
    var_lookup = load_isolate_variant_lookup(payload_dir / 'tables')
    # no longer needed
    # click.echo('Refactoring tables...')
    # refactor_sbj_history_core(payload_dir / 'tables', var_lookup)

    click.echo('Refactoring excludes...')
    refactor_sbj_history_core(payload_dir / 'excluded', var_lookup)
