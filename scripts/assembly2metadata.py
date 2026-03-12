#!/usr/bin/env python3
"""
assembly2metadata.py - Fetch BioSample metadata for NCBI Assembly accessions.

Given a list of NCBI Assembly accessions (GCA_* or GCF_*), this script resolves
each accession to its linked BioSample record via Entrez esummary, then fetches
the full BioSample metadata and writes a tab- or comma-delimited table.

Usage:
    python scripts/assembly2metadata.py -e your@email.com -i assemblies.txt -o output.tsv
    python scripts/assembly2metadata.py -e your@email.com -s GCA_000001405.15 GCF_000001635.27
    python scripts/assembly2metadata.py -e your@email.com -s GCA_000001405.15 -o output.csv

Key flags:
    -e / --email       Required. Email address for NCBI Entrez queries.
    -i / --input       Input file with one assembly accession per line.
                       Omit to read from stdin.
    -s / --assemblies  One or more assembly accessions on the command line.
    -o / --output      Output file (.tsv or .csv). Omit to write to stdout.
    -u / --update      Re-query accessions already present in the output file.
    --debug            Print verbose debug information to stderr.

Output:
    TSV (default) or CSV table with assembly_accession as the first column,
    followed by biosample_accession and all BioSample attribute fields found
    across the queried records. Missing fields are left blank.

If the output file already exists, existing records are loaded first and only
new accessions are queried (unless --update is set).

Requires: biopython
Install with: pip install biopython

NCBI rate limits: ~0.34 s sleep between requests (3 req/s unauthenticated).
Set Entrez.api_key for 10 req/s with a registered API key.
"""

import sys
import os
import re
import time
import csv
import argparse
from Bio import Entrez


parser = argparse.ArgumentParser(description='Extract BioSample metadata from NCBI Assembly accessions.')

parser.add_argument('-o', '--output', metavar=("out"), required=False,
                    help='Output table presenting lookup results'
                    ' if file already exists will update the file'
                    ' or append with new samples')

parser.add_argument('-i', '--input', metavar=("in"), required=False,
                    help='Input file of assembly accessions')

parser.add_argument('-s', '--assemblies', metavar=('assemblies'),
                    required=False, nargs='+',
                    help='Assembly accessions to look up')

parser.add_argument('-e', '--email', required=True,
                    help='Input your email address for Entrez queries')

parser.add_argument('-u', '--update', required=False, action='store_true',
                    help='Re-query samples already present in output file')

parser.add_argument('--debug', required=False, action='store_true',
                    help='Debug this running')


def fetch_assembly_info(assembly_acc):
    """Fetch assembly information including biosample accession."""
    try:
        handle = Entrez.esearch(db="assembly", term=assembly_acc, retmax=1)
        record = Entrez.read(handle)
        handle.close()

        if not record['IdList']:
            print(f"Warning: No assembly found for {assembly_acc}", file=sys.stderr)
            return None

        assembly_id = record['IdList'][0]

        handle = Entrez.esummary(db="assembly", id=assembly_id)
        summary = Entrez.read(handle, validate=False)
        handle.close()

        if 'DocumentSummarySet' in summary and 'DocumentSummary' in summary['DocumentSummarySet']:
            doc_summary = summary['DocumentSummarySet']['DocumentSummary'][0]
            biosample_acc = doc_summary.get('BioSampleAccn', '')
            return biosample_acc

        return None
    except Exception as e:
        print(f"Error fetching assembly {assembly_acc}: {e}", file=sys.stderr)
        return None


def fetch_biosample_metadata(biosample_acc):
    """Fetch all metadata fields from a biosample."""
    try:
        handle = Entrez.esearch(db="biosample", term=biosample_acc, retmax=1)
        record = Entrez.read(handle)
        handle.close()

        if not record['IdList']:
            print(f"Warning: No biosample found for {biosample_acc}", file=sys.stderr)
            return {}

        biosample_id = record['IdList'][0]

        handle = Entrez.efetch(db="biosample", id=biosample_id, rettype="xml")
        data = handle.read()
        handle.close()

        from xml.etree import ElementTree as ET
        root = ET.fromstring(data)

        metadata = {}
        for attr in root.findall('.//Attribute'):
            attr_name = attr.get('harmonized_name') or attr.get('attribute_name', 'unknown')
            metadata[attr_name] = attr.text if attr.text else ''

        return metadata
    except Exception as e:
        print(f"Error fetching biosample {biosample_acc}: {e}", file=sys.stderr)
        return {}


def main():
    args = parser.parse_args()
    Entrez.email = args.email

    separator = "\t"
    assemblies = {}
    header_set = set()

    if args.output:
        if re.search(r'\.(csv|CSV)$', args.output):
            separator = ","
        if os.path.exists(args.output):
            with open(args.output, "rt") as fh:
                csvin = csv.reader(fh, delimiter=separator)
                header = next(csvin)
                for h in header[1:]:
                    header_set.add(h)
                for row in csvin:
                    sampleid = row[0]
                    assemblies[sampleid] = {}
                    for i in range(1, len(header)):
                        assemblies[sampleid][header[i]] = row[i]

    query = set()
    if args.assemblies:
        for r in args.assemblies:
            if args.update or r not in assemblies:
                query.add(r)
    else:
        with open(args.input, 'rt') if args.input else sys.stdin as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if args.update or line not in assemblies:
                    if args.debug:
                        print(f"adding line {line}")
                    query.add(line)

    print(f"Processing {len(query)} assembly accessions...", file=sys.stderr)

    all_fields = set(['assembly_accession', 'biosample_accession']) | header_set

    for i, assembly_acc in enumerate(sorted(query), 1):
        print(f"Processing {i}/{len(query)}: {assembly_acc}", file=sys.stderr)

        biosample_acc = fetch_assembly_info(assembly_acc)
        time.sleep(0.34)  # NCBI rate limit: 3 requests per second

        if not biosample_acc:
            assemblies[assembly_acc] = {'biosample_accession': ''}
            continue

        metadata = fetch_biosample_metadata(biosample_acc)
        time.sleep(0.34)  # NCBI rate limit

        assemblies[assembly_acc] = {
            'biosample_accession': biosample_acc,
            **metadata
        }
        all_fields.update(metadata.keys())

    sorted_fields = sorted(all_fields - {'assembly_accession'})
    outheader = ['assembly_accession'] + sorted_fields

    with open(args.output, "wt") if args.output else sys.stdout as outfh:
        writer = csv.writer(outfh, delimiter=separator)
        writer.writerow(outheader)
        for acc in sorted(assemblies):
            row = [acc] + [assemblies[acc].get(f, '') for f in sorted_fields]
            writer.writerow(row)

    print(f"Total records: {len(assemblies)}", file=sys.stderr)
    print(f"Total metadata fields: {len(all_fields)}", file=sys.stderr)


if __name__ == "__main__":
    main()
