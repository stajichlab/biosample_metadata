# biosample_metadata
This tool supports extracting BioSample metadata into a user friendly table. It assumes you will provide a list of valid BioSamples and will write out a table.

If an existing table is provided it will be used to update that table with additional columns as present in the collection of samples provided or added to.

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

# Author(s)

Jason Stajich - jason.stajich[at]ucr.edu, http://lab.stajich.org
