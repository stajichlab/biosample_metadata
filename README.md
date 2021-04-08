# biosample_metadata
This tool supports extracting BioSample metadata into a user friendly table. It assumes you will provide a list of valid BioSamples and will write out a table. 

If an existing table is provided it will be used to update that table with additional columns as present in the collection of samples provided or added to. 

# Usage

```
biosample2table [ --in sample_list.txt] [--sample sampleid sampleid2 ...] --out result_table.tsv [--update]
``` 


| Argument | Description |
| --------- | -------- | 
| -o/--out |  Write out the table to this file, if the file already exists, parse it in and reuse the existing data, and add to it. Will use file extension (.tab/.tsv or .csv) to determine output format (default is tab delimited if extension doesn't match).  If `--out` is omitted will write to STDOUT (but will be unable to use previously cached results of course) |
| --update | use the existing table for sample ids, but also re-lookup the fields |
| -i/--in | Input file of biosamples, one per line for the query, if `-` is provided as filename will read in list from STDIN |
| -s/--sample | List of samples to query instead of providing as input file |
| --sra | Instead of BioSamples input IDs are SRR numbers |

# Author(s)
Jason Stajich - jason.stajich[at]ucr.edu, http://lab.stajich.org
