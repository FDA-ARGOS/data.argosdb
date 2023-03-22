import json
import click
import requests
from lxml import etree
from io import StringIO
from datetime import datetime
from pathlib import Path
from ..cli import cli
from ..utils.csvv import load_csv, dump_csv

INSERTIONS_URL = 'http://cov-glue-viz.cvr.gla.ac.uk/ins.php'
MUTATIONS_URL = 'http://cov-glue-viz.cvr.gla.ac.uk/mutations.php'
FIELDS = [
    'ref_name',
    'gene',
    'position',
    'amino_acid',
    'count',
    'proportion',
    'date_updated'
]


def convert_gene(gene):
    return {
        'ORF1ab/nsp1': 'nsp1',
        'ORF1ab/nsp2': 'nsp2',
        'ORF1ab/nsp3': 'PLpro',
        'ORF1ab/nsp4': 'nsp4',
        'ORF1ab/nsp5A-B': '_3CLpro',
        'ORF1ab/nsp6': 'nsp6',
        'ORF1ab/nsp7': 'nsp7',
        'ORF1ab/nsp8': 'nsp8',
        'ORF1ab/nsp9': 'nsp9',
        'ORF1ab/nsp10': 'nsp10',
        'ORF1ab/nsp12': 'RdRP',
        'ORF1ab/nsp13': 'nsp13',
        'ORF1ab/nsp14A2-B': 'nsp14',
        'ORF1ab/nsp15A1-B': 'nsp15',
        'ORF1ab/nsp16': 'nsp16',
        'S': 'S',
        'ORF3a': 'ORF3a',
        'E': 'E',
        'M': 'M',
        'ORF6': 'ORF6',
        'ORF7a': 'ORF7a',
        'ORF7b': 'ORF7b',
        'ORF8': 'ORF8',
        'N': 'N',
        'ORF10': 'ORF10'
    }.get(gene)


def convert_mut_row(row):
    gene = convert_gene(row['ORF'])
    if not gene:
        return
    pos = int(row['CodonNum'], 10)
    aa = row['Mutation'][-1]
    count = int(row['Count'], 10)
    proportion = float(row['Proportion'])
    date_updated = row['DateUpdated']
    return {
        'ref_name': 'Singer20',
        'gene': gene,
        'position': pos,
        'amino_acid': aa,
        'count': count,
        'proportion': proportion,
        'date_updated': date_updated
    }


def get_text(node):
    text = node.text or ''
    children = node.getchildren()
    for child in children:
        text += get_text(child)
    return text


def parse_response(resp, table_xpath):
    parser = etree.HTMLParser()
    tree = etree.parse(StringIO(resp.text), parser)
    updated_prefix = 'Updated: '
    updated, = tree.xpath(
        '//div[@class="navbar"]/a[starts-with(text(), {})][1]'
        .format(json.dumps(updated_prefix))
    )
    updated = datetime.strptime(
        updated.text[len(updated_prefix):],
        '%d/%m/%Y'
    ).date().isoformat()
    table, = tree.xpath(table_xpath)
    header = [h.text for h in table.xpath('./thead/tr/td')]
    rows = table.xpath('./tbody/tr')
    for row in rows:
        row = [get_text(d) for d in row.xpath('./td')]
        row = dict(zip(header, row))
        row['DateUpdated'] = updated
        yield row


def fetch_mutations(
    lineage_list='All',
    orf_list='All',
    mut_type_list='nonsyn',
    min_count_list=1,
    min_prop_list=0.01
):
    resp = requests.post(
        MUTATIONS_URL,
        data={
            'lineageList': lineage_list,
            'orfList': orf_list,
            'mutTypeList': mut_type_list,
            'minCountList': min_count_list,
            'minPropList': min_prop_list,
            'refreshTable': 'Refresh Table'
        }
    )
    rows = {}
    for row in parse_response(resp, '//table[@id="mutTable"][1]'):
        row = convert_mut_row(row)
        key = (row['gene'], row['position'], row['amino_acid'])
        rows.setdefault(key, {
            'ref_name': row['ref_name'],
            'gene': row['gene'],
            'position': row['position'],
            'amino_acid': row['amino_acid'],
            'count': 0,
            'proportion': .0,
            'date_updated': row['date_updated']
        })
        rows[key]['count'] += row['count']
        rows[key]['proportion'] += row['proportion']
    for row in rows.values():
        row['proportion'] = '{:f}'.format(row['proportion'])
        yield row


@cli.command()
@click.argument(
    'payload_dir',
    type=click.Path(
        dir_okay=True,
        exists=True,
        file_okay=False
    )
)
def update_glue_prevalence(payload_dir):
    payload_dir = Path(payload_dir)
    prev_dir = payload_dir / 'tables' / 'amino_acid_prevalence'
    prev_dir.mkdir(exist_ok=True)
    glue_csv = prev_dir / 'glue.csv'
    genes = load_csv(payload_dir / 'tables' / 'genes.csv')
    genes = {g['gene']: int(g['gene_order']) for g in genes}
    rows = sorted(
        fetch_mutations(min_prop_list='1E-05'),
        key=lambda r: (
            genes[r['gene']],
            r['position'],
            r['amino_acid']
        )
    )
    dump_csv(
        glue_csv,
        rows,
        FIELDS
    )
