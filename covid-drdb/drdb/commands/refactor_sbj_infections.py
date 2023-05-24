import click
from pathlib import Path
from ..cli import cli
from ..utils.csvv import load_csv, dump_csv, load_multiple_csvs


def load_immune_status(tables_dir):
    subjects = load_csv(tables_dir / 'subjects.csv')
    return {
        (sbj['ref_name'], sbj['subject_name']): sbj['immune_status']
        for sbj in subjects
    }


def load_severity(tables_dir):
    all_severity = load_multiple_csvs(tables_dir / 'subject_severity')
    lookup = {}
    for sbjsev in all_severity:
        sbjkey = (sbjsev['ref_name'], sbjsev['subject_name'])
        if sbjkey not in lookup:
            lookup[sbjkey] = []
        lookup[sbjkey].append(sbjsev)
    return lookup


def update_infections(tables_dir):
    # immune_status_lookup = load_immune_status(tables_dir)
    severity_lookup = load_severity(tables_dir)
    inf_dir = tables_dir / 'subject_infections'

    for inf_path in inf_dir.iterdir():
        sbjinfs = load_csv(inf_path)
        for sbjinf in sbjinfs:
            sbjkey = (sbjinf['ref_name'], sbjinf['subject_name'])
            # sbjinf['immune_status'] = immune_status_lookup[sbjkey]
            infdate = sbjinf['infection_date']
            for sev in severity_lookup.get(sbjkey, []):
                if infdate >= sev['start_date'] and infdate <= sev['end_date']:
                    sbjinf['severity'] = sev['severity']
                else:
                    print('unmatched: {}'.format(sbjkey))
        dump_csv(
            inf_path,
            records=sbjinfs,
            headers=[
                'ref_name',
                'subject_name',
                'infection_date_cmp',
                'infection_date',
                'infected_var_name',
                'location',
                'immune_status',
                'severity',
                'section'
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
def refactor_sbj_infections(payload_dir):
    payload_dir = Path(payload_dir)
    update_infections(payload_dir / 'excluded')
