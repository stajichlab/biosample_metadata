# biosample_metadata

Tools for retrieving BioSample and Assembly metadata from NCBI Entrez.

## Scripts

### `scripts/biosample2table.py`

Fetches metadata for BioSample or SRA accessions from NCBI and writes a TSV/CSV table.

```
python scripts/biosample2table.py -e your@email.com -i samples.txt -o output.tsv
python scripts/biosample2table.py -e your@email.com -s SAMN12345678 SAMN87654321
python scripts/biosample2table.py -e your@email.com --sra -s SRR1234567
```

**Key flags:**
- `-e` / `--email` — required; NCBI Entrez email
- `-i` / `--input` — input file (one accession per line); omit to read from stdin
- `-s` / `--sample` — space-separated accessions on the command line
- `-o` / `--output` — output file (`.tsv` or `.csv`); omit to write to stdout
- `-u` / `--update` — re-query accessions already present in the output file
- `--sra` — treat input accessions as SRR/ERR/DRR IDs (converts to BioSample first)
- `--showsra` — include SRA run accessions in output
- `--debug` — verbose debug output

If the output file already exists, existing records are loaded first; only new accessions are queried (unless `--update` is set).

---

### `scripts/assembly2metadata.py`

Fetches BioSample metadata for NCBI Assembly accessions (GCA_/GCF_) and writes a TSV/CSV table.

```
python scripts/assembly2metadata.py -e your@email.com -i assemblies.txt -o output.tsv
python scripts/assembly2metadata.py -e your@email.com -s GCA_000001405.15 GCF_000001635.27
```

**Key flags:**
- `-e` / `--email` — required; NCBI Entrez email
- `-i` / `--input` — input file (one assembly accession per line); omit to read from stdin
- `-s` / `--assemblies` — space-separated assembly accessions on the command line
- `-o` / `--output` — output file (`.tsv` or `.csv`); omit to write to stdout
- `-u` / `--update` — re-query accessions already present in the output file
- `--debug` — verbose debug output

Each assembly accession is resolved to a BioSample accession via `esummary`, then full metadata is fetched from the `biosample` database.

## Requirements

Install dependencies with:

```
pip install biopython
```

Or using the provided requirements file:

```
pip install -r requirements.txt
```

## NCBI Rate Limits

Both scripts sleep ~0.34 s between Entrez requests to respect the NCBI limit of 3 unauthenticated requests/second. Register an API key at https://www.ncbi.nlm.nih.gov/account/ and set `Entrez.api_key` to raise this to 10 requests/second.
