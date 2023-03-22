import click
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union
from collections import namedtuple
from ..cli import cli
from ..utils.csvv import load_csv, dump_csv
from ..utils.questionary_option import QuestionaryOption

optional = namedtuple('optional', 'tpl_name')


def load_articles(payload_dir: Path) -> List[Dict]:
    articles_csv: Path = payload_dir / 'tables' / 'articles.csv'
    return load_csv(articles_csv)


def abort_if_ref_name_used(
    ctx: click.Context,
    param: click.Option,
    value: str
) -> str:
    articles: List[Dict] = load_articles(Path(ctx.params['payload_dir']))
    for a in articles:
        if a['ref_name'] == value:
            raise click.BadParameter(
                'Specified {} {} is already used by article {}'
                .format(param.opts[-1], a['ref_name'], a['doi'] or a['url']),
                ctx,
                param
            )
    return value


def abort_if_ref_name_unused(
    ctx: click.Context,
    param: click.Option,
    value: str
) -> str:
    articles: List[Dict] = load_articles(Path(ctx.params['payload_dir']))
    if not any(a['ref_name'] == value for a in articles):
        raise click.BadParameter(
            'Specified {} {} does not exist'
            .format(param.opts[-1], value),
            ctx,
            param
        )
    return value


def add_article(
    ref_name: str,
    doi: str,
    first_author: str,
    year: str,
    payload_dir: Path
) -> None:
    articles: List[Dict] = load_articles(payload_dir)
    articles.append({
        'ref_name': ref_name,
        'doi': doi if doi.startswith('10.') else None,
        'url': None if doi.startswith('10.') else doi,
        'first_author': first_author,
        'year': year,
        'date_added': datetime.now().strftime('%Y-%m-%d')
    })
    dump_csv(
        payload_dir / 'tables' / 'articles.csv',
        records=articles,
        headers=list(articles[0].keys()),
        BOM=True
    )


def save_xx_tpl(
    ref_name: str,
    tpl_dir: Path,
    filename_suffix: str,
    repeat: int,
    **extra_fields: Optional[str]
) -> Path:
    tpl_path: Path = tpl_dir / '{}{}.csv'.format(
        ref_name.lower(), filename_suffix)
    ex_path: Path = next(
        tpl_dir.glob('*{}.csv'.format(filename_suffix))
    )
    ex_rows: List[Dict] = load_csv(ex_path)
    dump_csv(
        tpl_path,
        records=[
            {**extra_fields, 'ref_name': ref_name}
        ] * repeat,
        headers=list(ex_rows[0].keys()),
        BOM=True
    )
    return tpl_path


@cli.command()
@click.argument(
    'payload_dir',
    is_eager=True,
    type=click.Path(
        dir_okay=True,
        exists=True,
        file_okay=False
    )
)
@click.option(
    '--ref-name', prompt="Enter refName",
    callback=abort_if_ref_name_used,
    help='Reference name of article to be entered')
@click.option(
    '--doi', prompt="Enter DOI/URL",
    help='DOI of article to be entered')
@click.option(
    '--first_author', prompt="Enter first author (e.g. Harari, S)",
    help='DOI of article to be entered')
@click.option(
    '--year', prompt="Enter publication year (e.g. 2022)",
    help='Publication year of article to be entered')
@click.option(
    '--study-type', prompt='Select study type',
    type=click.Choice([
        'plasma-susceptibility',
        'mab-susceptibility',
        'drug-susceptibility',
        'susceptibility',
        'invivo-selection',
        'invitro-selection'
    ], case_sensitive=False),
    cls=QuestionaryOption,
    help='Study type'
)
def new_study(
    payload_dir: str,
    ref_name: str,
    doi: str,
    first_author: str,
    year: str,
    study_type: str
) -> None:
    tables_dir: Path = Path(payload_dir) / 'tables'

    tpls: Dict[str, Tuple[str, str, int, Dict[str, Optional[str]]]] = {
        'isolates': (
            'isolates.d',
            '-iso',
            4,
            {
                'iso_name': 'hCoV-19/...',
                'site_directed': 'FALSE',
                'gisaid_id': '1234567',
                'expandable': 'TRUE'
            }
        ),
        'isolate_mutations': (
            'isolate_mutations.d',
            '-isomuts',
            10,
            {'iso_name': 'hCoV-19/...', 'gene': 'S'}
        ),
        'subject_isolates': (
            'subject_isolates',
            '-sbjiso',
            4,
            {
                'iso_name': 'hCoV-19/...',
                'subject_name': 'Patient ...',
                'iso_source': 'NP',
                'iso_culture': 'FALSE',
                'collection_date_cmp': '='
            }
        ),
        'subject_infections': (
            'subject_infections',
            '-inf',
            2,
            {
                'subject_name': 'Patient ...',
                'infection_date_cmp': '<',
            }
        ),
        'subject_treatments': (
            'subject_treatments',
            '-prx',
            4,
            {
                'subject_name': 'Patient ...',
                'start_date_cmp': '<',
                'end_date_cmp': '>'
            }
        ),
        'subject_plasma': (
            'subject_plasma',
            '-plasma',
            4,
            {
                'subject_name': 'Patient ...',
                'collection_date_cmp': '~'
            }
        ),
        'rx_antibodies': (
            'rx_antibodies',
            '-mab',
            2,
            {}
        ),
        'rx_compounds': (
            'rx_compounds',
            '-drug',
            2,
            {}
        ),
        'rx_potency': (
            'rx_potency',
            '-pot',
            8,
            {
                'iso_name': 'hCoV-19/...',
                'date_added': datetime.now().strftime('%Y-%m-%d')
            }
        ),
        'rx_fold': (
            'rx_fold',
            '-fold',
            4,
            {
                'control_iso_name': 'B.1 full genome',
                'iso_name': 'hCoV-19/...',
                'date_added': datetime.now().strftime('%Y-%m-%d')
            }
        ),
        'ref_isolate_pairs': (
            'ref_isolate_pairs',
            '-pair',
            2,
            {
                'control_iso_name': 'B.1 full genome',
                'iso_name': 'hCoV-19/...',
            }
        ),
        'subject_vaccines': (
            'subject_vaccines',
            '-vacc',
            2,
            {
                'subject_name': 'Patient ...',
                'vaccination_date_cmp': '~'
            }
        )
    }

    study_type_tpls: Dict[str, List[Union[str, optional]]] = {
        'susceptibility': [
            'isolates',
            'isolate_mutations',
            'ref_isolate_pairs',
            optional('subject_infections'),
            optional('subject_vaccines'),
            optional('subject_plasma'),
            optional('rx_antibodies'),
            optional('rx_compounds'),
            optional('rx_potency'),
            optional('rx_fold'),
        ],
        'plasma-susceptibility': [
            'isolates',
            'isolate_mutations',
            'ref_isolate_pairs',
            optional('subject_infections'),
            optional('subject_vaccines'),
            'subject_plasma',
            optional('rx_potency'),
            optional('rx_fold'),
        ],
        'mab-susceptibility': [
            'isolates',
            'isolate_mutations',
            'ref_isolate_pairs',
            'rx_antibodies',
            optional('rx_potency'),
            optional('rx_fold'),
        ],
        'drug-susceptibility': [
            'isolates',
            'isolate_mutations',
            'ref_isolate_pairs',
            'rx_compounds',
            optional('rx_potency'),
            optional('rx_fold'),
        ],
        'invivo-selection': [
            'isolates',
            'isolate_mutations',
            'subject_isolates',
            'subject_infections',
            optional('subject_plasma'),
            optional('rx_antibodies'),
            optional('rx_compounds'),
            'subject_treatments',
            optional('subject_vaccines'),
        ],
        'invitro-selection': [
            'invitro_selection_results'
        ]
    }

    click.echo('=' * 50)
    click.echo('  Study type: {}'.format(study_type))
    click.echo('=' * 50)
    click.echo()
    click.echo(
        'Following template files are created under `{}` directory. '
        'You should edit & update them for entering this new study. '
        'You can delete the optional files if not used.'
        .format(payload_dir)
    )

    add_article(ref_name, doi, first_author, year, Path(payload_dir))

    tpl_dir: str
    fname_suffix: str
    repeat: int
    extra: Dict[str, Optional[str]]
    for tplkey in study_type_tpls[study_type]:
        opt: bool = False
        if isinstance(tplkey, optional):
            opt = True
            tpl_name = tplkey.tpl_name
        else:
            tpl_name = tplkey
        tpl_dir, fname_suffix, repeat, extra = tpls[tpl_name]
        tpl = save_xx_tpl(
            ref_name,
            tables_dir / tpl_dir,
            fname_suffix,
            repeat,
            **extra
        )
        click.echo('- {}{}'.format(tpl, ' (optional)' if opt else ''))
