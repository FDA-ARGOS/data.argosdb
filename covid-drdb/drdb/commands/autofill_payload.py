import re
import click
from pathlib import Path
from more_itertools import unique_everseen
from operator import itemgetter
from typing import List, Dict, Optional, Tuple

from ..cli import cli
from ..utils.gene_position import translate_gene_position
from ..utils.csvv import (
    load_csv,
    dump_csv,
    load_multiple_csvs,
    CSVReaderRow,
    CSVWriterRow
)


def autofill_invitros(tables_dir: Path) -> None:
    rows: List[CSVReaderRow]
    invitro: Path
    invitros: Path = tables_dir / 'invitro_selection_results'
    for invitro in invitros.iterdir():
        if invitro.suffix.lower() != '.csv':
            click.echo('Skip {}'.format(invitro))
            continue
        rows = load_csv(invitro)
        click.echo('Write to {}'.format(invitro))
        dump_csv(
            invitro,
            records=rows,
            headers=[
                'ref_name',
                'rx_name',
                'backbone',
                'gene',
                'position',
                'amino_acid',
                'section',
                'date_added'
            ],
            BOM=True
        )


def autofill_invivos(tables_dir: Path) -> None:
    row: CSVReaderRow
    rows: List[CSVReaderRow]
    invivo: Path
    invivos: Path = tables_dir / 'ref_invivo'
    for invivo in invivos.iterdir():
        if invivo.suffix.lower() != '.csv':
            click.echo('Skip {}'.format(invivo))
            continue
        rows = load_csv(invivo)

        for row in rows:
            if not row.get('note'):
                row['note'] = None

        click.echo('Write to {}'.format(invivo))
        dump_csv(
            invivo,
            records=rows,
            headers=[
                'ref_name',
                'subject_name',
                'collection_date',
                'section',
                'note',
                'date_added',
            ],
            BOM=True
        )


def autofill_rx(tables_dir: Path) -> None:
    file_path: Path

    rxmabs: List[CSVReaderRow] = load_multiple_csvs(
        tables_dir / 'rx_antibodies')

    rxps: List[
        Dict[str, Optional[str]]
    ] = load_multiple_csvs(tables_dir / 'subject_plasma')

    rxdrugs: List[CSVReaderRow] = load_multiple_csvs(
        tables_dir / 'rx_compounds')

    unclassified_rx: List[CSVReaderRow] = []
    file_path = tables_dir / 'unclassified-rx.csv'
    if file_path.exists():
        unclassified_rx = load_csv(file_path)

    # warning: do not add rx_names from invitro, invivo, or DMS tables Only the
    # rx_antibodies, subject_plasma and unclassified_rx should be the source of
    # treatments table.

    treatments: List[CSVWriterRow] = list(unique_everseen([
        {'ref_name': rx['ref_name'],
         'rx_name': rx['rx_name']}
        for rx in rxmabs + rxps + rxdrugs + unclassified_rx
    ]))
    click.echo('Write to {}'.format(tables_dir / 'treatments.csv'))
    dump_csv(
        tables_dir / 'treatments.csv',
        records=treatments,
        headers=['ref_name', 'rx_name'],
        BOM=True
    )


def autofill_sbj_plasma(tables_dir: Path) -> None:
    rows: List[CSVReaderRow]
    rxcp: Path
    rxcps: Path = tables_dir / 'subject_plasma'
    for rxcp in rxcps.iterdir():
        if rxcp.suffix.lower() != '.csv':
            click.echo('Skip {}'.format(rxcp))
            continue
        rows = load_csv(rxcp)
        for row in rows:
            if not row.get('collection_date_cmp'):
                row['collection_date_cmp'] = '='
            if not row.get('subject_name'):
                row['subject_name'] = row['rx_name']
            if not row.get('location'):
                row['location'] = None
            if not row.get('section'):
                row['section'] = None
        click.echo('Write to {}'.format(rxcp))
        dump_csv(
            rxcp,
            records=rows,
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


def autofill_dms(tables_dir: Path) -> None:
    row: CSVReaderRow
    ace2_binding: Path = tables_dir / 'dms' / 'dms_ace2_binding.csv'
    ace2_contact: Path = tables_dir / 'dms' / 'dms_ace2_contact.csv'
    rows: List[CSVReaderRow] = load_csv(ace2_binding)

    pos_ace2_contact = [
        int(rec['position'])
        for rec in load_csv(ace2_contact)
        if rec['ace2_contact']
    ]

    for row in rows:
        if not row['ace2_binding']:
            row['ace2_binding'] = None
        if not row['expression']:
            row['expression'] = None
        row['ace2_contact'] = (
            'TRUE' if int(row['position']) in pos_ace2_contact else 'FALSE')

    click.echo('Write to {}'.format(ace2_binding))
    dump_csv(
        ace2_binding,
        records=rows,
        headers=[
            # 'ref_name',
            'gene',
            'position',
            'amino_acid',
            'ace2_binding',
            'expression',
            'ace2_contact',
            # 'ace2_assay',
            # 'backbone',
        ],
        BOM=True,
    )

    escape_score: Path = tables_dir / 'dms' / 'dms_escape_results.csv'
    rows = load_csv(escape_score)

    for row in rows:
        if not row['escape_score']:
            row['escape_score'] = None

    click.echo('Write to {}'.format(escape_score))
    dump_csv(
        escape_score,
        records=rows,
        headers=[
            'ref_name',
            'rx_name',
            'gene',
            'position',
            'amino_acid',
            'escape_score',
            'method',
            # 'backbone',
        ],
        BOM=True,
    )


def sort_csv(file_path: Path, *key_list: str) -> None:
    records: List[CSVReaderRow] = load_csv(file_path)
    records.sort(key=itemgetter(*key_list))
    dump_csv(file_path, records)


def autofill_subjects(tables_dir: Path) -> None:
    known_subjects: Dict[Tuple[str, str], CSVReaderRow] = {
        (r['ref_name'], r['subject_name']): r
        for r in load_csv(tables_dir / 'subjects.csv')
        if r['ref_name'] is not None and
        r['subject_name'] is not None
    }

    sbj_plasma: List[CSVReaderRow] = load_multiple_csvs(
        tables_dir / 'subject_plasma')
    sbj_isolates: List[CSVReaderRow] = load_multiple_csvs(
        tables_dir / 'subject_isolates')

    subjects: List[CSVWriterRow] = sorted(
        unique_everseen([
            {'ref_name': rx['ref_name'],
             'subject_name': rx['subject_name'],
             'subject_species': (
                 known_subjects
                 .get((rx['ref_name'], rx['subject_name']), {})
                 .get('subject_species') or 'Human'
             ),
             'birth_year': (
                 known_subjects
                 .get((rx['ref_name'], rx['subject_name']), {})
                 .get('birth_year') or 'NULL'
             ),
             'num_subjects': (
                 known_subjects
                 .get((rx['ref_name'], rx['subject_name']), {})
                 .get('num_subjects') or 1
             ),
             }
            for rx in sbj_plasma + sbj_isolates
            if rx['ref_name'] is not None and
            rx['subject_name'] is not None

        ]),
        key=lambda rx: (rx['ref_name'] or '', rx['subject_name'] or '')
    )
    click.echo('Write to {}'.format(tables_dir / 'subjects.csv'))
    dump_csv(
        tables_dir / 'subjects.csv',
        records=subjects,
        headers=[
            'ref_name',
            'subject_name',
            'subject_species',
            'birth_year',
            'num_subjects'
        ],
        BOM=True
    )


def autofill_sbj_treatments(tables_dir: Path) -> None:
    row: CSVReaderRow
    rows: List[CSVReaderRow]
    prx: Path
    prx_list: Path = tables_dir / 'subject_treatments'
    for prx in prx_list.iterdir():
        if prx.suffix.lower() != '.csv':
            click.echo('Skip {}'.format(prx))
            continue
        rows = load_csv(prx)
        for row in rows:
            if not row.get('start_date_cmp'):
                row['start_date_cmp'] = '='
            if not row.get('end_date_cmp'):
                row['end_date_cmp'] = '='
            if not row.get('section'):
                row['section'] = None
        click.echo('Write to {}'.format(prx))
        dump_csv(
            prx,
            records=rows,
            headers=[
                'ref_name',
                'subject_name',
                'rx_name',
                'start_date_cmp',
                'start_date',
                'end_date_cmp',
                'end_date',
                'dosage',
                'dosage_unit',
                'section'
            ],
            BOM=True
        )


def autofill_sbj_infections(tables_dir: Path) -> None:
    row: CSVReaderRow
    rows: List[CSVReaderRow]
    pth: Path
    pth_list: Path = tables_dir / 'subject_infections'
    for pth in pth_list.iterdir():
        if pth.suffix.lower() != '.csv':
            click.echo('Skip {}'.format(pth))
            continue
        rows = load_csv(pth)
        for row in rows:
            if not row.get('infection_date_cmp'):
                row['infection_date_cmp'] = '='
            if not row.get('infected_var_name'):
                row['infected_var_name'] = 'Unknown variant'
            if not row.get('location'):
                row['location'] = None
            if not row.get('section'):
                row['section'] = None
            if not row.get('immune_status'):
                row['immune_status'] = None
            if not row.get('severity'):
                row['severity'] = None
        click.echo('Write to {}'.format(pth))
        dump_csv(
            pth,
            records=rows,
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


def autofill_isolates(tables_dir: Path) -> None:
    row: CSVReaderRow
    rows: List[CSVReaderRow]
    iso: Path
    iso_list: Path = tables_dir / 'isolates.d'
    for iso in [tables_dir / 'isolates.csv', *iso_list.iterdir()]:
        if iso.suffix.lower() != '.csv':
            click.echo('Skip {}'.format(iso))
            continue
        rows = load_csv(iso)
        for row in rows:
            if not row.get('var_name'):
                row['var_name'] = None
            if row.get('site_directed') != 'TRUE':
                row['site_directed'] = 'FALSE'
            if not row.get('gisaid_id'):
                row['gisaid_id'] = None
            if not row.get('genbank_accn'):
                row['genbank_accn'] = None
            if not row.get('sra_accn'):
                row['sra_accn'] = None
            if row.get('expandable') != 'FALSE':
                row['expandable'] = 'TRUE'

            gisaid = row.get('gisaid_id')
            if gisaid and gisaid.startswith('EPI_ISL_'):
                row['gisaid_id'] = gisaid[8:]
            genbank = row.get('genbank_accn')
            sra = row.get('sra_accn')
            if genbank and re.search(r'^[CES]RR|^SAMN', genbank):
                sra = genbank
                genbank = None
            row['genbank_accn'] = genbank
            row['sra_accn'] = sra
        dump_csv(
            iso,
            records=rows,
            headers=[
                'iso_name',
                'var_name',
                'site_directed',
                'gisaid_id',
                'genbank_accn',
                'sra_accn',
                'expandable'
            ],
            BOM=True
        )


def autofill_isomuts(tables_dir: Path) -> None:
    row: CSVReaderRow
    rows: List[CSVReaderRow]
    isomuts: Path
    isomuts_list: Path = tables_dir / 'isolate_mutations.d'
    for isomuts in [
        tables_dir / 'isolate_mutations.csv',
        *isomuts_list.iterdir()
    ]:
        if isomuts.suffix.lower() != '.csv':
            click.echo('Skip {}'.format(isomuts))
            continue
        rows = load_csv(isomuts)
        for row in rows:
            gene: Optional[str] = row['gene']
            pos_text: Optional[str] = row['position']
            if not gene or not pos_text:
                continue
            gene, pos = translate_gene_position(gene, int(pos_text))
            aa: Optional[str] = row['amino_acid']
            if not aa:
                continue
            if aa == '-':
                aa = 'del'
            elif aa == '_' or '_' in aa:
                aa = 'ins'
            elif aa == '*':
                aa = 'stop'
            row['gene'] = gene
            row['position'] = str(pos)
            row['amino_acid'] = aa
        dump_csv(
            isomuts,
            records=rows,
            headers=[
                'iso_name',
                'gene',
                'position',
                'amino_acid',
                'count',
                'total'
            ],
            BOM=True
        )


def autofill_sbj_isolates(tables_dir: Path) -> None:
    row: CSVReaderRow
    rows: List[CSVReaderRow]
    pth: Path
    pth_list: Path = tables_dir / 'subject_isolates'
    for pth in pth_list.iterdir():
        if pth.suffix.lower() != '.csv':
            click.echo('Skip {}'.format(pth))
            continue
        rows = load_csv(pth)
        for row in rows:
            if not row.get('collection_date_cmp'):
                row['collection_date_cmp'] = '='
            if not row.get('iso_source'):
                row['iso_source'] = 'Pending'
            if not row.get('iso_culture'):
                row['iso_culture'] = 'FALSE'
            if not row.get('location'):
                row['location'] = None
            if not row.get('section'):
                row['section'] = None
        click.echo('Write to {}'.format(pth))
        dump_csv(
            pth,
            records=rows,
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


def autofill_sbj_vaccines(tables_dir: Path) -> None:
    row: CSVReaderRow
    rows: List[CSVReaderRow]
    pth: Path
    pth_list: Path = tables_dir / 'subject_vaccines'
    for pth in pth_list.iterdir():
        if pth.suffix.lower() != '.csv':
            click.echo('Skip {}'.format(pth))
            continue
        rows = load_csv(pth)
        for row in rows:
            if not row.get('vaccination_date_cmp'):
                row['vaccination_date_cmp'] = '='
            if not row.get('location'):
                row['location'] = None
            if not row.get('section'):
                row['section'] = None
        click.echo('Write to {}'.format(pth))
        dump_csv(
            pth,
            records=rows,
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


def autofill_rx_potency(tables_dir: Path) -> None:
    row: CSVReaderRow
    rows: List[CSVReaderRow]
    pth: Path
    pth_list: Path = tables_dir / 'rx_potency'
    for pth in pth_list.iterdir():
        if pth.suffix.lower() != '.csv':
            click.echo('Skip {}'.format(pth))
            continue
        rows = load_csv(pth)
        for row in rows:
            row['potency_type'] = row.get('potency_type') or 'NT50'
            row['cumulative_count'] = row.get('cumulative_count') or '1'
        click.echo('Write to {}'.format(pth))
        dump_csv(
            pth,
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


def autofill_rx_fold(tables_dir: Path) -> None:
    row: CSVReaderRow
    rows: List[CSVReaderRow]
    susc: Path
    suscs: Path = tables_dir / 'rx_fold'
    for susc in suscs.iterdir():
        if susc.suffix.lower() != '.csv':
            click.echo('Skip {}'.format(susc))
            continue
        rows = load_csv(susc)
        for row in rows:
            if not row.get('fold'):
                row['fold'] = None
                row['fold_cmp'] = None
            elif not row.get('fold_cmp'):
                row['fold_cmp'] = '='
            if not row.get('cumulative_count'):
                row['cumulative_count'] = '1'
            if not row.get('resistance_level'):
                row['resistance_level'] = None
            if not row.get('assay_name'):
                row['assay_name'] = None
            if not row.get('control_iso_name'):
                row['control_iso_name'] = 'Control'
            if not row.get('ineffective'):
                row['ineffective'] = None
            if not row.get('potency_type'):
                row['potency_type'] = 'IC50'
            if row.get('potency_type') == '0':
                row['potency_type'] = 'IC50'
            if row.get('potency_type') == '90':
                row['potency_type'] = 'IC90'
            if row.get('potency_type') == '50':
                row['potency_type'] = 'IC50'

        click.echo('Write to {}'.format(susc))
        dump_csv(
            susc,
            records=rows,
            headers=[
                'ref_name',
                'rx_name',
                'control_iso_name',
                'iso_name',
                'section',
                'assay_name',
                'potency_type',
                'fold_cmp',
                'fold',
                'resistance_level',
                'ineffective',
                'cumulative_count',
                'date_added'
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
def autofill_payload(payload_dir: str) -> None:
    payload_dir_path = Path(payload_dir)
    tables_dir: Path = payload_dir_path / 'tables'
    antibodies: Path = tables_dir / 'antibodies.csv'
    antibody_targets: Path = tables_dir / 'antibody_targets.csv'
    autofill_rx(tables_dir)
    autofill_invitros(tables_dir)
    autofill_invivos(tables_dir)
    autofill_sbj_plasma(tables_dir)
    autofill_dms(tables_dir)
    autofill_rx_fold(tables_dir)
    autofill_rx_potency(tables_dir)
    autofill_isolates(tables_dir)
    autofill_isomuts(tables_dir)

    autofill_subjects(tables_dir)
    autofill_sbj_infections(tables_dir)
    autofill_sbj_isolates(tables_dir)
    autofill_sbj_vaccines(tables_dir)
    autofill_sbj_treatments(tables_dir)

    sort_csv(antibodies, 'ab_name')
    sort_csv(antibody_targets, 'ab_name')
