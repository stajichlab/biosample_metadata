# Changelog

## 2026-03-12

### `scripts/assembly2metadata.py` — major fixes

The script had several bugs that prevented it from running at all:

- **Added missing imports**: `argparse`, `re`, `os` were used but not imported, causing immediate `NameError` crashes.
- **Added `--assemblies` / `-s` argument**: `args.assemblies` was referenced but the argument was never defined in the parser.
- **Rewrote `main()`**: The original `main()` ignored all argparse setup — it re-read `sys.argv[1]` directly, never set `Entrez.email`, and hardcoded the output filename. It now uses `args` throughout.
- **Honored `--output`**: Output file was hardcoded as `assembly_biosample_metadata.tsv`; now uses `args.output` (defaults to stdout).
- **Removed duplicate `time.sleep`**: Each loop iteration slept twice (0.34 s × 2 = 0.68 s); reduced to one sleep per Entrez call.
- **Removed unused `defaultdict` import**.
- **Prefer `harmonized_name` attribute**: BioSample XML attributes now use `harmonized_name` when available (consistent with `biosample2table.py`).
- **Skip blank lines** in input file.

### `scripts/biosample2table.py` — minor fixes

- **Removed unused `pprint` import** and the unused `pp` object.
- **Guarded "adding line" print** with `if args.debug:` — it was printing unconditionally to stdout on every input line.
