import re
import csv
import time
import json
import click
import zipfile
import requests
from io import BytesIO
from pathlib import Path
from typing import Dict, Any, List, Tuple

from tqdm import tqdm  # type: ignore

from ..cli import cli
from ..utils.csvv import dump_csv

DOWNLOADER_PATTERN = re.compile(
    r'(/downloader\.php\?file_name=epitope_table_export_\d+\.zip)'
)
FILENAME_PATTERN = re.compile(r'(epitope_table_export_\d+)')
COMPILED_DATA_PATTERN = re.compile(r'var compiledData\s*=\s*({[^\n]+})\s*;\n')


def query(
    params: Dict[str, Any],
    worker_type: int,
    encrypted: str
) -> Any:
    handler: str = 'http://www.iedb.org/WorkRequestHandler.php'
    result: Dict[str, Any]
    resp: requests.Response = requests.post(
        handler,
        data={
            'worker_type': worker_type,
            'encrypted': encrypted,
            'params': json.dumps(params)
        }
    )
    params = resp.json()
    while worker_type == 4:
        time.sleep(1)
        resp = requests.get(handler, params=params)
        result = resp.json()
        if result['result_code'] == 'SUCCESS':
            return json.loads(result['result_data'])

    result_url: List[str]
    while worker_type == 2:
        time.sleep(1)
        resp = requests.get(
            'http://www.iedb.org/exportDownload/{id}'.format(**params))
        result_url = DOWNLOADER_PATTERN.findall(resp.text)
        if result_url:
            resp = requests.get('http://www.iedb.org' + result_url[0])
            fp = BytesIO(resp.content)
            zipfp: zipfile.ZipFile = zipfile.ZipFile(fp)
            filename = FILENAME_PATTERN.findall(result_url[0])
            if filename:
                csvfp = zipfp.open(filename[0] + '.csv')
                csvtext = csvfp.read().decode('UTF-8').splitlines()
                header = [
                    h1 + ' ' + h2
                    for h1, h2 in zip(*list(csv.reader(csvtext[:2])))
                ]
                return list(csv.DictReader([','.join(header)] + csvtext[2:]))


def fetch_csv() -> Tuple[List[str], List[Dict[str, Any]]]:
    result: Dict[str, Any]
    params: Dict[str, Any] = {
        'filter_view': 'epitope',
        'structure_type': 'linear_sequence',
        'e_value': 'exact',
        'source_organism_ids[]': [
            'http://purl.obolibrary.org/obo/NCBITaxon_2697049'
        ],
        'source_antigen_ids[]': ['90584'],
        'pos_measured_response': 1,
        'assay_t_cell_all': 1,
        'assay_b_cell_all': 0,
        'assay_mhc_all': 0,
        'mhc_class_type': 'class_I',
        'mhc_class_ids[]': ['-1'],
        'mhc_resolution': 'any',
        'host_organism_type': 'human_host',
        'host_organism_ids[]': [
            'http://purl.obolibrary.org/obo/NCBITaxon_9606'
        ],
        'disease_type': 'any_disease',
        'reference_type': 'anyreference',
        'receptor_type': 'any_type',
        'receptor_chain_type': 'any_type',
        'receptor_chain_sequence_type': 'cdr3',
        'receptor_chain_e_value': 'exact',
        'tcr_receptor_type': 'any_type',
        'tcr_receptor_chain_type': 'any_type',
        'tcr_receptor_chain_sequence_type': 'cdr3',
        'tcr_receptor_chain_e_value': 'exact',
        'tcell_only_assay_type': 'any_assay_type',
        'bcell_only_assay_type': 'any_assay_type',
        'ligand_only_assay_type': 'any_assay_type',
        'bcr_receptor_type': 'any_type',
        'bcr_receptor_chain_type': 'any_type',
        'bcr_receptor_chain_sequence_type': 'cdr3',
        'bcr_receptor_chain_e_value': 'exact',
        'count': '',
        'page_num': '',
        'sort_col': 'ref_count',
        'sort_dir': 'desc',
        'items_per_page': '1000',
        'start_char': '',
        'sort_col_type': 'num',
        'list_type': 'peptidic',
        'search_type': 'simple_search'
    }
    params_export: Dict[str, Any] = {
        **params,
        'export_type': 'base_export',
        'bread_crumb_values': {}
    }
    result = query(params, worker_type=4, encrypted='false')
    params_export['result_count'] = result['count']
    header: List[str] = [h['data_key'] for h in result['list_header']]
    rows: List[Dict[str, Any]] = [
        {key: val['value'] for key, val in row.items()}
        for row in result['results']
    ]

    addrows = query(params_export, worker_type=2, encrypted='true')
    addrows_lookup = {
        row['Epitope Epitope ID']: {
            'structure_id': row['Epitope Epitope ID'],
            'structure_description_sort': row['Epitope Description'],
            'structure_start': row['Epitope Starting Position'],
            'structure_end': row['Epitope Ending Position']
        }
        for row in addrows
    }
    for row in tqdm(rows):
        struct_id = row['structure_id']
        addrow = addrows_lookup[struct_id]
        if addrow['structure_description_sort'] != \
                row['structure_description_sort']:
            click.echo('Epitope are not the same: ' +
                       row['structure_description_sort'] + ' and ' +
                       addrow['structure_description_sort'], err=True)
            raise click.Abort()
        row.update(addrow)
        resp = requests.get(
            'http://www.iedb.org/epitope/{}'.format(struct_id)
        )
        compiled_data = COMPILED_DATA_PATTERN.findall(resp.text)
        mhc_molecule = []
        if compiled_data:
            data = json.loads(compiled_data[0])
            for tbl in data['data']:
                for item in tbl['data']:
                    if 'mhc_molecule' in item:
                        mhc_molecule.append(item['mhc_molecule'])
        row['mhc_molecule'] = '; '.join(mhc_molecule)
    return header + ['structure_start', 'structure_end', 'mhc_molecule'], rows


@cli.command()
@click.argument(
    'payload_dir',
    type=click.Path(
        dir_okay=True,
        exists=True,
        file_okay=False
    )
)
def fetch_iedb_epitopes(payload_dir: str) -> None:
    payload_dir_path = Path(payload_dir)
    iedb_epitopes: Path = (
        payload_dir_path / 'excluded' / 'iedb_t_cell_epitopes.csv'
    )
    header, rows = fetch_csv()
    dump_csv(iedb_epitopes, rows, header)
