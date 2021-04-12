#!/usr/bin/env python3
import os, re, csv, sys
import argparse

from Bio import Entrez
import xml.etree.ElementTree as ET
import pprint

pp = pprint.PrettyPrinter(indent=4)

def indent(elem, level=0):
    i = "\n" + level*"  "
    j = "\n" + (level-1)*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for subelem in elem:
            indent(subelem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = j
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = j
    return elem


parser = argparse.ArgumentParser(description='Extract BioSample metadata from NCBI Entrez.')

parser.add_argument('-o', '--output', metavar=("out"), required = False,
                    help='Output table presenting lookup results, if file already exists will update the file or append with new samples')
parser.add_argument('-i', '--input', metavar=("in"),required = False,
                    help='Input file of sample names')

parser.add_argument('-s', '--sample', metavar=('samples'), required = False, nargs='+',
                    help='Input file of sample names')

parser.add_argument('-e', '--email', required = True,
                    help='Input your email address for Entrez queries')

parser.add_argument('-u', '--update', required = False, action='store_true',
                    help='Input sample names as list')

parser.add_argument('--sra', required = False, action='store_true',
                    help='Sample names are SRR IDs not BioSamples')


parser.add_argument('--debug', required = False, action='store_true',
                    help='Debug this running')


args = parser.parse_args()

Entrez.email = args.email

biosamples = {}

separator = "\t"
query = set()
header_set = set()
seenSRR = set()
if args.output:
    if re.search(r'\.(csv|CSV)$',args.output):
        separator = ","
    if os.path.exists(args.output):
        with open(args.output,"rt") as fh:
            csvin = csv.reader(fh,delimiter=separator)
            header = next(csvin)
            for h in header[1:]:
                header_set.add(h)
            for row in csvin:
                sampleid = row[0]
                biosamples[sampleid] = {}
                # skip column 1 since it is the ID column
                for i in range(1,len(header)):
                    biosamples[sampleid][header[i]] = row[i]
                    if header[i] == "SRA_Run":
                        seenSRR.add(row[i])

#  input list of sample ids can be either
# from -s samplelist
# or by file (--in)
# or by stdin
if args.sample:
    for r in args.sample:
        if args.update or not (r in biosamples or r in seenSRR):
            query.add(r)
else:
    with open(args.input, 'rt') if args.input else sys.stdin as fh:
        for line in fh:
            line = line.strip()
            if args.update or not (line in biosamples or line in seenSRR):
                print("adding line {}".format(line))
                query.add(line)


sampid2sra = {}

# if the sra parameter is set then we need to convert these IDs to BioSample IDs
if args.sra:
    newquery = set()
    for sraid in query:
        if sraid.startswith("SAM"):
            print("skipping {} as it looks like a biosample, are you sure you meant to use --sra".format(sraid))
            continue
        handle = Entrez.efetch(db="sra", id=sraid)
        tree = ET.parse(handle)
        root = tree.getroot()

        for sd in root.iter('SAMPLE'):
            if args.debug:
                indent(sd)
                ET.dump(sd)
            if "accession" in sd.attrib:
                if args.debug:
                    print("SRA BioSample Accession is {}".format(sd.attrib["accession"]))
                newquery.add(sd.attrib["accession"])
                if args.debug:
                    print("{} => {}".format(sd.attrib["accession"],sraid))
                sampid2sra[sd.attrib["accession"]] = sraid
                sampid2sra[ sraid ] = sd.attrib["accession"]
    query = newquery # update the query set with this

# we need to convert SAMNXXX to internal IDs for biosamples that NCBI uses
# requires using esearch
# this loop will finish with a unique list of samples to run final lookup on
# in the sampidquery set
sampidquery = set()
biosamp2sra = {}
for qname in query:
    if args.debug:
        print("query is {}".format(qname))
    m = re.match(r'^(SAM|[SED]RS)',qname)
    if not m:
        print("query is {} and looks like it isn't a biosample, did you forget to specify --sra".format(qname))
        continue

    handle = Entrez.esearch(db="biosample", term=qname)
    tree = ET.parse(handle)
    root = tree.getroot()
    if args.debug:
        indent(root)
        ET.dump(root)

    for idlist in root.iter("IdList"):
        for id in idlist:
            if args.debug:
                print("ID found {}",format(id))
                indent(id)
                ET.dump(id)
            if id.tag == "Id":
                if qname in sampid2sra:
                    biosamp2sra[ id.text ] = sampid2sra[qname]
                sampidquery.add(id.text)
# now query internal NCBI sampid values with efetch
# and get back all the metadata columns
# which are associated with the attributes values
for sampid in sampidquery:
    handle = Entrez.efetch(db="biosample", id=sampid)
    if args.debug:
        print("sampid is {}".format(sampid))
    tree = ET.parse(handle)
    root = tree.getroot()
    for sample in root:
        if args.debug:
            indent(sample)
            ET.dump(sample)
        BIOSAMPLE = sample.attrib['accession']
        if BIOSAMPLE not in biosamples:
            biosamples[BIOSAMPLE] = {}
        if sampid in biosamp2sra:
            biosamples[BIOSAMPLE]['SRA_Run'] = biosamp2sra[sampid]
            biosamples[BIOSAMPLE]['SRA_SampID'] = sampid2sra [ biosamp2sra[sampid] ]
            header_set.add('SRA_Run')
            header_set.add('SRA_SampID')
        for attributes in root.iter('Attributes'): # sample to root?
            for metadata in attributes:
                keyname = metadata.attrib['attribute_name']
                if 'harmonized_name' in metadata.attrib:
                    keyname = metadata.attrib['harmonized_name']
                header_set.add(keyname)
                biosamples[BIOSAMPLE][keyname] = metadata.text

with open(args.output,"wt") if args.output else sys.stdout as outfh:
    outcsv = csv.writer(outfh, delimiter=separator)
    outheader = ['BioSample']
    sorted_header_set = sorted(header_set)
    outheader.extend(sorted_header_set)
    outcsv.writerow(outheader)
    for sample in sorted(biosamples):
        outrow = [sample]
        for h in sorted_header_set:
            if h in biosamples[sample]:
                outrow.append(biosamples[sample][h])
            else:
                outrow.append('')
        outcsv.writerow(outrow)
