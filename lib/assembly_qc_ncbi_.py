#!/usr/bin/env python3
"""Assembly QC NCBI

This script will parse an XML for a single feature and return a list via terminal
or file, if an output is supplied.
"""
__version__ = "0.1.0"
__status__ = "BETA"

import os
import sys
import csv
from ftplib import FTP
import xml.etree.ElementTree as ET
from argparse import ArgumentParser, SUPPRESS
from Bio import Entrez


def usr_args():
    """User Arguments

    User supplied arguments from command line for function

    Returns
    -------
        ArgumentParser objects to be digested by subsequent functions.
    """

    parser = ArgumentParser(
        add_help=False,
        prog='data_sheet_validator',
        description="Release Note Generator. "
            "Used to generate release notes from two differnt versions.")

    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument('-f', '--file',
        required=True,
        help="File and path for parsing.")

    # required.add_argument('-', '--new',
    #     required=True,
    #     help="Directory for new version of dictionary")

    optional.add_argument('-o', '--output',
        help="Output file. If no output is provided, the resulting table will \
            be output to the terminal.")
    optional.add_argument('-v', '--version',
        action='version',
        version='%(prog)s ' + __version__)
    optional.add_argument('-h', '--help',
        action='help',
        default=SUPPRESS,
        help='show this help message and exit')

    if len(sys.argv) <= 1:
        sys.argv.append('--help')

    return parser.parse_args()

def read_datasheet(input_path):
    """Read input data sheet

    Designed to take the ARGOS ngs ID List and Selection Criteria and extract
    all unique NCBI assembly accessions

    Parameters
    ----------
    input_path: str
        file path/name to be parsed
    """

    genome_assembly_ids =[]

    with open(input_path, 'r', encoding='utf') as accessions:
        data = csv.reader(accessions, delimiter='\t')
        next(data)
        for row in data:
            if row[2] not in genome_assembly_ids and row[2] != '':
                genome_assembly_ids.append(row[2])

    genomes = ','.join(genome_assembly_ids)
    efetch = f'efetch -db assembly -id {genomes} -format docsum > home/genomes.xml'
    os.system(efetch)

def parse_xml(xml_file, samples):
    """Parse XML file

    Parameters
    ----------
    xml_file: str
        file path/name to be parsed
    output: str, optional
        file path/name for data to be output to
    """
    genome_stats = {}
    genome_assembly_id = lineage = taxonomy_id = bco_id = schema_version \
        = analysis_platform = analysis_platform_object_id = bioproject \
        = biosample = sra_run_id = ngs_read_file_name = ngs_read_file_source \
        = organism_name = ref_genome_acc = genomic_section \
        = num_chromosomes = num_genes = num_segments  = assembly_gc_content \
        = length = n50  = n75 = n90 = l50  = l75 = number_of_n \
        = percent_assembly_greater_20x = query_coverage_against_reference \
        = percent_identity_against_reference = reads_unaligned = id_method \
        = assembly_type = assembly_level = assembly_score = '-'

    analysis_platform = 'NCBI'
    bco_id = 'https://biocomputeobject.org/ARGOS_000038/DRAFT'
    schema_version = 'v0.9'
    genomic_section = 'all'

    root = ET.parse(xml_file).getroot()
    for item in root:
        if item.tag == "DbBuild":
            continue
        for summary in item.findall('./'):
            if summary.tag == 'AssemblyAccession':
                genome_assembly_id = summary.text
            if summary.tag == 'Id':
                analysis_platform_object_id = summary.text
            if summary.tag == 'Organism':
                organism_name = summary.text
            if summary.tag == 'Taxid':
                taxonomy_id = summary.text
                lineage = get_lineage(taxonomy_id)
            if summary.tag == 'AssemblyType':
                assembly_type = summary.text
            if summary.tag == 'GB_BioProjects':
                if summary[0][0].tag == 'BioprojectAccn':
                    bioproject = summary[0][0].tag
            if summary.tag == 'RS_BioProjects':
                if summary[0][0].tag == 'BioprojectAccn':
                    bioproject = summary[0][0].text
            if summary.tag == 'BioSampleAccn':
                biosample = summary.text
            if summary.tag == 'FtpPath_GenBank':
                ngs_read_file_source = summary.text + '_genomic.gbff.gz'
                ngs_read_file_name = summary.text.split('/')[-1] \
                    + '_genomic.gbff.gz'
            if summary.tag == 'Synonym':
                for synonm in summary:
                    if synonm.tag == 'RefSeq':
                        ref_genome_acc = synonm.text
            if summary.tag == 'Meta':
                for meta in summary.findall('./'):
                    if meta.tag == 'Stats':
                        for stat in meta.findall('./'):
                            if stat.attrib['category'] == 'chromosome_count':
                                num_chromosomes = stat.text
                            if stat.attrib['category'] == 'total_length':
                                length = stat.text
                            if stat.attrib['category'] == 'contig_n50':
                                n50 = stat.text
                            if stat.attrib['category'] == 'contig_l50':
                                l50 = stat.text
                    if meta.tag == 'assembly-status':
                        assembly_level = meta.text
            if summary.tag == 'FtpPath_Stats_rpt':
                print(summary.text)
                genome_stats = stats_report(summary.text, genome_stats)
                num_segments = genome_stats['top-level-count']
                assembly_gc_content = genome_stats['gc-perc']
                n75 = genome_stats['scaffold-N75']
                n90 = genome_stats['scaffold-N90']

        samples[genome_assembly_id] = lineage , taxonomy_id , bco_id , schema_version \
            , analysis_platform , analysis_platform_object_id , bioproject \
            , biosample , sra_run_id , ngs_read_file_name , ngs_read_file_source \
            , organism_name , ref_genome_acc , genomic_section \
            , num_chromosomes , num_genes , num_segments  , assembly_gc_content \
            , length , n50  , n75 , n90 , l50  , l75 , number_of_n \
            , percent_assembly_greater_20x , query_coverage_against_reference \
            , percent_identity_against_reference , reads_unaligned , id_method \
            , assembly_type , assembly_level , assembly_score

    return samples

def stats_report(stats_ftp, genome_stats):
    """Assembly Stats Report

    """

    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    stats_ftp = stats_ftp.split('://')[1]
    stats_file = stats_ftp.split('/')[-1]
    genome_dir = stats_ftp.split(stats_file)[0]
    genome_dir = genome_dir.replace('ftp.ncbi.nlm.nih.gov/', '')
    ftp.cwd(genome_dir)
    with open('home/' + stats_file, 'wb') as file:
        ftp.retrbinary(f'RETR {stats_file}', file.write )

    with open('home/' + stats_file, 'r', encoding='utf8') as statfile:
        genome_tsv = csv.reader(statfile, delimiter="\t")
        for row in genome_tsv:
            if row[0].startswith('#'):
                continue
            if row[0] == 'all' and row[4] == 'gc-perc':
                genome_stats['gc-perc'] = row[5]
            if row[0] == 'all' and row[4] == 'scaffold-count':
                genome_stats['top-level-count'] = row[5]
            if row[0] == 'all' and row[4] == 'scaffold-N75':
                genome_stats['scaffold-N75'] = row[5]
            if row[0] == 'all' and row[4] == 'scaffold-N90':
                genome_stats['scaffold-N90'] = row[5]

    return genome_stats

def get_lineage(taxonomy_id):
    """Get lineage

    Retrieves the NCBI Taxonomic lineage given a tax id
    """
    Entrez.email = "hadley_king@gwu.edu"
    Entrez.api_key = os.getenv('NCBI_API_KEY')
    handle = Entrez.efetch(db="taxonomy", id=taxonomy_id)
    record = Entrez.read(handle)
    return record[0]['Lineage']

def sample_output( samples, header, output):
    """Sample Output

    If an output file is supplied in the user arguments the samples dictionary
    will be output in TSV to the supplied file. If no output is provided, the
    samples dictionary will print to the terminal in a TSV format.

    Parameters
    ----------
    header: lst of str
        List of column headers for the output
    samples: dict
        dictionary of samples
    outputs: str, optional
        file path/name to output data to
    """

    if output:
        sample_file = os.path.abspath(output)
        with open(sample_file, 'w',  encoding='utf8') as file:
            writer = csv.writer(file, header, delimiter='\t')
            writer.writerow(header)
            for key in samples:
                row = [key]
                for item in range(len(samples[key])):
                    row.append(samples[key][item])
                writer.writerow(row)
        print('Samples written to ', sample_file)
    else:
        print('\t'.join(item for item in header))
        for run in samples:
            print(run, '\t', '\t'.join(str(item) for item in samples[run]))

def main():
    """Main"""
    samples = {}
    header = ['genome_assembly_id', 'lineage', 'taxonomy_id', 'bco_id', 'schema_version', \
        'analysis_platform', 'analysis_platform_object_id', 'bioproject', \
        'biosample',  'sra_run_id', 'ngs_read_file_name', 'ngs_read_file_source', \
        'organism_name', 'ref_genome_acc', 'genomic_section', \
        'num_chromosomes', 'num_genes', 'num_segments', 'assembly_gc_content', \
        'length', 'n50', 'n75', 'n90', 'l50', 'l75', 'number_of_n',  \
        'percent_assembly_greater_20x', 'query_coverage_against_reference', \
        'percent_identity_against_reference',  'reads_unaligned', 'id_method',\
        'assembly_type', 'assembly_level', 'assembly_score']

    args = usr_args()
    read_datasheet(args.file)
    parse_xml('home/genomes.xml', samples)
    sample_output(samples, header, args.output)

if __name__ == "__main__":
    main()
