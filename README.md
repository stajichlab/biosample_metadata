# Get Metadata for Biosamples
This tool supports extracting BioSample metadata into a user friendly table. It assumes you will provide a list of valid BioSamples and will write out a table.

If an existing table is provided it will be used to update that table with additional columns as present in the collection of samples provided or added to.

# Requirements

* [Biopython](http://biopython.org)

# Usage

```
biosample2table.py [ --in sample_list.txt] [--sample sampleid sampleid2 ...] --out result_table.tsv [--update]
```

## Retrieve a single biosample

Lookup a single Biosample and output tab delimited file
```
./scripts/biosample2table.py -s SAMN12327137 --out result_table.tsv -e yourname@gmail.com
```

Results will be
```
BioSample	collection_date	geo_loc_name	host	isolate	sample_type	strain
SAMN12327137	2006	USA:Sierra Nevada Mtns	Rana sierrae	Bd_JAM81	isolate	GPL
```

## Output comma delimited file
```
./scripts/biosample2table.py -s SAMN12327137 --out result_table.csv -e yourname@gmail.com
```
BioSample,collection_date,geo_loc_name,host,isolate,sample_type,strain
SAMN12327137,2006,USA:Sierra Nevada Mtns,Rana sierrae,Bd_JAM81,isolate,GPL

## Lookup a single SRA record and find its BioSample
```
./scripts/biosample2table.py -s SRR14174621 --sra --out result_table.tsv -e yourname@gmail.com
```

Results will be
```
BioSample,SRA_Run,SRA_SampID,Sample No.,age,biomaterial_provider,isolate,sex,tissue,treatment
SAMN18650164,SRR14174621,SRS8658905,17,54,Liverpool Hospital,Diabetic Patient,male,Skin,Midpoint
```

## Lookup list of IDs from a single filename

If you have a list of IDs either SRA or BioSample you can use the `--in` option

```
./scripts/biosample2table.py --in samplefile.txt --out result_table.tsv -e yourname@gmail.com
```

# Using STDIN and STDOUT

If you had a list of ids you wanted to pass from another program uyou can pass that in and omit the `--in` option

```
echo SAMN18650164 | ./scripts/biosample2table.py --out result_table.tsv -e yourname@gmail.com
```

Omitting the `--out` option will print the results out to STDOUT

```
echo SAMN18650164 | ./scripts/biosample2table.py -e yourname@gmail.com
```


Mixing `--in` (or STDIN with no `--in`) and the `-s` will prefer sample list input provided by the `-s` option and ignore any stdin or `--in` file input.

```
echo SAMN18650164 | ./scripts/biosample2table.py -e yourname@gmail.com
```

# Command line arguments

| Argument | Description |
| --------- | -------- |
| -o/--out |  Write out the table to this file, if the file already exists, parse it in and reuse the existing data, and add to it. Will use file extension (.tab/.tsv or .csv) to determine output format (default is tab delimited if extension doesn't match).  If `--out` is omitted will write to STDOUT (but will be unable to use previously cached results of course) |
| --update | use the existing table for sample ids, but also re-lookup the fields |
| -i/--in | Input file of biosamples, one per line for the query, if `-` is provided as filename will read in list from STDIN |
| -s/--sample | List of samples to query instead of providing as input file |
| --sra | Instead of BioSamples input IDs are SRR numbers |

# Assembly Accession to Metadata (`assembly2metadata.py`)

`scripts/assembly2metadata.py` fetches BioSample metadata for NCBI Assembly accessions
(GCA_* or GCF_*). Each assembly is resolved to its linked BioSample via Entrez, then the
full metadata is retrieved and written to a table.

## Usage

```
assembly2metadata.py -e your@email.com [-i assemblies.txt] [-s GCA_... GCF_...] [-o output.tsv] [--update]
```

## Retrieve metadata for a single assembly

```
./scripts/assembly2metadata.py -s GCA_000001405.15 -o result_table.tsv -e yourname@gmail.com
```

Results will be a tab-delimited table with `assembly_accession` as the first column,
`biosample_accession` as the second, followed by all BioSample attribute fields:

```
assembly_accession	biosample_accession	collection_date	geo_loc_name	...
GCA_000001405.15	SAMN02981236	...	...	...
```

## Retrieve metadata for multiple assemblies

```
./scripts/assembly2metadata.py -s GCA_000001405.15 GCF_000001635.27 -o result_table.tsv -e yourname@gmail.com
```

## Retrieve from a file list

One assembly accession per line:

```
./scripts/assembly2metadata.py -i assemblies.txt -o result_table.tsv -e yourname@gmail.com
```

## Output comma-delimited file

Use a `.csv` extension to switch to comma-delimited output:

```
./scripts/assembly2metadata.py -s GCA_000001405.15 -o result_table.csv -e yourname@gmail.com
```

## Update existing results

Re-query assemblies already present in the output file:

```
./scripts/assembly2metadata.py -i assemblies.txt -o result_table.tsv -e yourname@gmail.com --update
```

## Using STDIN and STDOUT

```
echo GCA_000001405.15 | ./scripts/assembly2metadata.py -e yourname@gmail.com
```

Omitting `-o` writes results to stdout.

## Command line arguments

| Argument | Description |
| --------- | -------- |
| -e/--email | Required. Email address for NCBI Entrez queries. |
| -o/--output | Output file (.tsv or .csv). If the file already exists, existing records are reused and new ones are appended. Omit to write to stdout. |
| -i/--input | Input file of assembly accessions, one per line. Omit to read from stdin. |
| -s/--assemblies | One or more assembly accessions on the command line. |
| --update | Re-query accessions already present in the output file. |
| --debug | Print verbose debug information to stderr. |

# Author(s)

Jason Stajich - jason.stajich[at]ucr.edu, http://lab.stajich.org
