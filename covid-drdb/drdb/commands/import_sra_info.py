import click
from pathlib import Path
from lxml import etree
from Bio import Entrez  # type: ignore
from typing import Dict, List, Optional

from ..cli import cli
from ..utils.csvv import load_csv, dump_csv
from .new_study import abort_if_ref_name_unused

Entrez.email = 'philiptz@stanford.edu'

ISO_SOURCE_MAP: Dict[str, str] = {
    'nasopharyngeal': 'NP',
    'nasopharyngeal swab': 'NP',
    'endotracheal aspirate': 'ETA',
    'bronchoalveolar lavage': 'BAL',
    'oropharyngeal': 'OP',
    'oropharyngeal swab': 'OP',
    'plasma': 'Plasma',
    'saliva': 'Saliva',
    'stool': 'Stool'
}


def load_sra_info_lookup(
    accessions: List[str]
) -> Dict[str, Dict[str, Optional[str]]]:
    sra_lookup: Dict[str, Dict[str, Optional[str]]] = {}
    with Entrez.efetch(db='sra', id=accessions) as handle:
        root = etree.parse(handle)
        exps = root.xpath('/EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE')
        if not isinstance(exps, list):
            raise click.ClickException('Malformed XML received from Entrez')
        for exp in exps:
            if not isinstance(exp, etree._Element):
                continue
            samples = exp.xpath('./SAMPLE')
            if not isinstance(samples, list):
                continue
            sample = samples[0]
            if not isinstance(sample, etree._Element):
                continue
            attrs = sample.xpath('./SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE')
            if not isinstance(attrs, list):
                continue
            attr_lookup: Dict[str, Optional[str]] = {}
            for attr in attrs:
                if not isinstance(attr, etree._Element):
                    continue
                tags = attr.xpath('./TAG/text()')
                values = attr.xpath('./VALUE/text()')
                if not isinstance(tags, list) or not tags:
                    continue
                if not isinstance(values, list) or not values:
                    continue
                tag, value = tags[0], values[0]
                if isinstance(tag, str) and isinstance(value, str):
                    attr_lookup[tag] = value
            runs = exp.xpath('./RUN_SET/RUN/@accession')
            if not isinstance(runs, list):
                continue
            for run in runs:
                if not isinstance(run, str):
                    continue
                sra_lookup[run] = attr_lookup
    return sra_lookup


def update_all(
    isolates: List[Dict],
    isomuts: List[Dict],
    sbjisos: List[Dict],
    sra_lookup: Dict[str, Dict[str, Optional[str]]],
    update_iso_name: bool
) -> None:
    for isolate in isolates:
        accn: Optional[str] = isolate['genbank_accn']
        if not accn:
            continue
        attr_lookup: Optional[
            Dict[str, Optional[str]]
        ] = sra_lookup.get(accn)
        if not attr_lookup:
            continue
        orig_iso_name: str = isolate['iso_name']

        collection_date: Optional[str] = attr_lookup.get('collection_date')
        if collection_date:
            for item in sbjisos:
                if item['iso_name'] == orig_iso_name:
                    item['collection_date'] = collection_date

        iso_source: Optional[str] = attr_lookup.get('isolation_source')
        if iso_source:
            iso_source = ISO_SOURCE_MAP.get(iso_source.lower(), iso_source)
            for item in sbjisos:
                if item['iso_name'] == orig_iso_name:
                    item['iso_source'] = iso_source

        location: Optional[str] = attr_lookup.get('geo_loc_name')
        if location:
            for item in sbjisos:
                if item['iso_name'] == orig_iso_name:
                    item['location'] = location

        sra_iso_name: Optional[str] = attr_lookup.get('isolate')
        if update_iso_name and sra_iso_name:
            isolate['iso_name'] = sra_iso_name
            for item in isomuts + sbjisos:
                if item['iso_name'] == orig_iso_name:
                    item['iso_name'] = sra_iso_name


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
    '--ref-name', prompt="Enter new study's refName",
    callback=abort_if_ref_name_unused,
    help='Reference name of article to be entered')
@click.option(
    '--update-iso-name/--no-update-iso-name',
    is_flag=True,
    default=True,
    prompt='Should the `iso_name` be updated according to SRA?',
    help='Should the `iso_name` be updated according to SRA?'
)
def import_sra_info(
    payload_dir: str,
    ref_name: str,
    update_iso_name: bool
) -> None:
    tables_dir: Path = Path(payload_dir) / 'tables'
    lower_ref_name: str = ref_name.lower()
    iso_path: Path = (
        tables_dir / 'isolates.d' /
        '{}-iso.csv'.format(lower_ref_name)
    )
    isolates: List[Dict] = load_csv(iso_path)

    isomut_path: Path = (
        tables_dir / 'isolate_mutations.d' /
        '{}-isomuts.csv'.format(lower_ref_name)
    )
    isomuts: List[Dict] = load_csv(isomut_path)

    sbjiso_path: Path = (
        tables_dir / 'subject_isolates' /
        '{}-sbjiso.csv'.format(lower_ref_name)
    )
    sbjisos: List[Dict] = []
    if sbjiso_path.is_file():
        sbjisos = load_csv(sbjiso_path)

    accessions: List[str] = [
        iso['genbank_accn'] for iso in isolates
        if iso['genbank_accn']
    ]
    sra_lookup: Dict[
        str,
        Dict[str, Optional[str]]
    ] = load_sra_info_lookup(accessions)
    update_all(
        isolates,
        isomuts,
        sbjisos,
        sra_lookup,
        update_iso_name
    )
    dump_csv(
        iso_path,
        records=isolates,
        headers=list(isolates[0].keys()),
        BOM=True)
    click.echo('Write {}'.format(iso_path))
    dump_csv(
        isomut_path,
        records=isomuts,
        headers=list(isomuts[0].keys()),
        BOM=True)
    click.echo('Write {}'.format(isomut_path))
    dump_csv(
        sbjiso_path,
        records=sbjisos,
        headers=list(sbjisos[0].keys()),
        BOM=True)
    click.echo('Write {}'.format(sbjiso_path))
