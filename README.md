# Sjogren public omics preliminary aggregate results

This standalone snapshot contains complete locally computed aggregate statistics and pathway families from 17 accepted result packages. It includes age, culture, VST, captured-mixture, methylation and mechanistic reference contexts. Package names, methods and table roles preserve their original exploratory and reference scopes. This is not a complete Sjogren database, independent-cohort collection, clinical validation or causal conclusion. Scientific completion is false.

`manifest.json` lists the exact tables, counts, schemas, original/published hashes and retained identifiers. Current indexed pathways and historical/supporting pathways have separate status values. Every submitted row, missing-value literal, filter state, effect, P value and correction remains unchanged. Export only changes compression. No new fits, corrections, orthology or donor claims are introduced.

Python 3.9 or newer, standard library only. Run from any location:

```sh
python -B query.py status
python -B query.py tables --limit 10
python -B query.py gene CXCL13 --limit 10
python -B query.py gene PD-1 --limit 10
python -B query.py feature ENSG00000156234 --limit 10
python -B query.py pathway GO:0006955 --limit 10
python -B query.py validate
```

Queries stream compressed CSV and may scan millions of rows. Use `--package` or exact `--table` to narrow scope. Default output is 100 rows; `--offset` enables paging and `--all` explicitly streams all hits. Gene lookup is literal and casefolded with declared source multi-symbol delimiters; PD-1/PD-L1 map explicitly to PDCD1/CD274. Feature lookup is exact. Pathways match exact casefolded source IDs/names. Output gives table ID, original row number, query mapping and unmodified source fields. Missing hits do not establish measurement absence or a negative biological result. No numeric conversion or significance filtering occurs.

`methods.json` preserves contract scales, correction families, method evidence and accepted scope. Unknown tissue, n, confounding or exact pathway eligibility remain null/unknown; numeric suffixes are not inferred biological sample sizes. `source_registry.csv` preserves 200 source-status records as bounded context; these are not 200 independent cohorts or an exhaustive search. Historical query63 and its older comparison/target summaries are excluded. Raw inputs, individual metadata, matrices, barcode libraries and author-reported supplements are excluded.

Validation checks hashes, CSV headers, row widths and complete counts offline. Producer source equivalence and independent review are separately documented in `verification.json`. Repository visibility is controlled by its owner; these files do not declare a public release or scientific publication readiness.

One methylation logical table is stored in ordered, header-preserving gzip parts of at most 4 MiB for the observed API transport constraint. This changes physical packaging only: all 485,577 rows, literal values and original source row numbers remain unchanged. The manifest records contiguous ranges and the reconstructed original decompressed SHA256. Queries and validation treat these parts as one logical table. There remain 200 logical tables and 3,308,952 rows.
