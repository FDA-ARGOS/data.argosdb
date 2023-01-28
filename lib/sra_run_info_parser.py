#!/usr/bin/env python3
""" SRA Run Info Parser

This script will parse the XML result file from the SraRunInfo.xml. It takes
a directory as an input and will iterate through all files, parsing each one
to first a dictionary, and then will output to the terminal in
TSV format, or to a file if an output file is provided.
"""

import os
import sys
import csv
import argparse
import xml.etree.ElementTree as ET
__version__ = "0.1.0"
__status__ = "Beta"

def usr_args() -> argparse.ArgumentParser():
    """Program Arguments

    All arguments for process are defined here.
    """

    parser = argparse.ArgumentParser()

    parser = argparse.ArgumentParser(
        prog='sra_trace_xml.py',
        usage='%(prog)s [options]')


    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s ' + __version__)

    parser.add_argument('-d', '--directory',
        required=True,
        help="Directory for list of input files. The input files should be SRR \
            XML files.")
    parser.add_argument('-o', '--output',
        help="Output file. If no output is provided, the resulting table will be \
            output to the terminal.")
    if len(sys.argv) <= 1:
        sys.argv.append('--help')

    return parser.parse_args()

def parse_xml( xml_file: str, samples: dict)-> dict:
    """Parse XML file

    Parameters
    ----------
    xml_file: str
        file path/name to be parsed
    samples: dict
        dictionary of samples

    Returns
    -------
    samples: dict
        dictionary of samples
    """
    bco_id = 'ARGOS_000009'
    schema_version = 'v1.2'
    sra_run_id = xml_file.split('/')[-1].split('.')[0]

    organism_name = base_total = infraspecific_name = lineage \
        = genome_assembly_id = taxonomy_id = analysis_platform \
        = analysis_platform_object_id = bioproject = biosample = strain \
        = ngs_read_file_name = ngs_read_file_source \
        = ngs_gc_content = count_a = count_t = count_g \
        = count_c = count_n = instrument = id_method = wgs_accession \
        = strategy = n_content = ngs_score = ''

    avg_phred_score = max_read_length = min_read_length\
        = num_reads_unique = pos_outlier_count = codon_table = percent_coding\
        = percent_not_coding = density_n_per_read = complexity_percent\
        = non_complex_percent = avg_quality_a = avg_quality_t = avg_quality_g\
        = avg_quality_c = avg_read_length ='HIVE'

    try:
        root = ET.parse(xml_file).getroot()
        for item in root.findall('.//'):
            if item.tag == 'EXPERIMENT':
                analysis_platform_object_id = item.attrib['alias']
            if item.tag == 'STUDY':
                bioproject = item.attrib['alias']
        for exp in root.findall('.//EXPERIMENT/PLATFORM'):
            analysis_platform = exp[0].tag
            instrument = exp[0].tag
        for lib in root.findall('.//EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/'):
            if lib.tag == 'LIBRARY_STRATEGY':
                strategy = lib.text
        for run in root.findall('.//RUN'):
            for bases in run.findall('./Bases'):
                base_total = float(bases.attrib['count'])
                for base in bases.findall('./'):
                    if base.attrib['value'] == 'A':
                        count_a = float(base.attrib['count'])
                    if base.attrib['value'] == 'C':
                        count_c = float(base.attrib['count'])
                    if base.attrib['value'] == 'G':
                        count_g = float(base.attrib['count'])
                    if base.attrib['value'] == 'T':
                        count_t = float(base.attrib['count'])
                    if base.attrib['value'] == 'N':
                        count_n = float(base.attrib['count'])
            for item in run.findall('./SRAFiles/SRAFile'):
                if item.attrib['filename'] == run.attrib['accession']:
                    ngs_read_file_name = item.attrib['filename']
                    ngs_read_file_source = item.attrib['url']
            for item in run.findall('./IDENTIFIERS/PRIMARY_ID'):
                sra_run_id = item.text
        for sample in root.findall('.//SAMPLE'):
            strain = sample.attrib['alias']
            biosample = sample.attrib['accession']
        for sample in root.findall('.//SAMPLE/DESCRIPTION'):
            infraspecific_name = sample.text
        for name in root.findall('.//SAMPLE/SAMPLE_NAME/TAXON_ID'):
            taxonomy_id = name.text
        for name in root.findall('.//SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME'):
            organism_name = name.text
        for attribute in root.findall('.//SAMPLE/SAMPLE_ATTRIBUTES/'):
            if attribute[0].text in ('identification_method', 'identification method'):
                id_method = attribute[1].text
        ngs_gc_content = round((((count_g+count_c)/base_total)*100), 2)
        n_content = ((count_n/base_total)*100)

    except ET.ParseError:
        print(sra_run_id)
        organism_name = f'{sra_run_id}.xml had no data'

    samples[sra_run_id] = [organism_name, infraspecific_name, lineage,
        genome_assembly_id, taxonomy_id, bco_id, schema_version,
        analysis_platform, analysis_platform_object_id, bioproject,
        biosample, strain, ngs_read_file_name, ngs_read_file_source,
        ngs_gc_content, n_content, avg_phred_score, avg_read_length,
        max_read_length, min_read_length, num_reads_unique,
        pos_outlier_count, codon_table, percent_coding,
        percent_not_coding, density_n_per_read, complexity_percent,
        non_complex_percent, avg_quality_a, avg_quality_t,
        avg_quality_g, avg_quality_c, count_a, count_t, count_g,
        count_c, count_n, instrument, id_method, wgs_accession, strategy,
        ngs_score]

    return samples

def sample_output(samples, header, output):
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
        with open(sample_file, 'w', encoding='utf8') as file:
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
    """Main Function
    """
    samples = {}

    header = ['sra_run_id','organism_name', 'infraspecific_name', 'lineage',
    'genome_assembly_id', 'taxonomy_id', 'bco_id', 'schema_version',
    'analysis_platform', 'analysis_platform_object_id', 'bioproject',
    'biosample', 'strain', 'ngs_read_file_name', 'ngs_read_file_source',
    'ngs_gc_content', 'n_content', 'avg_phred_score', 'avg_read_length',
    'max_read_length', 'min_read_length', 'num_reads_unique',
    'pos_outlier_count', 'codon_table', 'percent_coding', 'percent_not_coding',
    'density_n_per_read', 'complexity_percent', 'non_complex_percent',
    'avg_quality_a', 'avg_quality_t', 'avg_quality_g', 'avg_quality_c',
    'count_a', 'count_t', 'count_g', 'count_c', 'count_n', 'instrument',
    'id_method', 'wgs_accession', 'strategy', 'ngs_score']

    args = usr_args()
    for item in os.listdir(args.directory):
        xml_file = os.path.join(args.directory, item)
        parse_xml(xml_file, samples)
    sample_output(samples, header, args.output)

if __name__ == "__main__":
    main()
