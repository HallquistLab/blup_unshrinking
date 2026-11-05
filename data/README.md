# REPAIR simulation data

This directory documents the Parquet files used by the REPAIR simulation
studies. The binary datasets are housed on Zenodo: together they are 
pproximately 3.6 GB. The files are *simulated data*, not human data.

## Access

```text
REPAIR simulation data DOI: 10.5281/zenodo.23965546
```

## Files

Note that due to the size of the files, download times may be longer.

| File | Download | Size (bytes) |
|---|---|---:|---|---|
| `study1v2.parquet` | [Zenodo](https://zenodo.org/records/23965547/files/study1v2.parquet) | 782478809 |
| `study2v2.parquet` | [Zenodo](https://zenodo.org/records/23965547/files/study2v2.parquet) | 984683117 |
| `study3v2.parquet` | [Zenodo](https://zenodo.org/records/23965547/files/study3v2.parquet) | 1920730132 |
| `study4v2.parquet` | [Zenodo](https://zenodo.org/records/23965547/files/study4v2.parquet) | 216399016 |

Verify a download with:

```sh
shasum -a 256 study1v2.parquet
shasum -a 256 study2v2.parquet
shasum -a 256 study3v2.parquet
shasum -a 256 study4v2.parquet
```

The files are Parquet tables and can be read with R's `arrow` package or
Python's `pyarrow` package:

```r
library(arrow)
study2 <- read_parquet("study2v2.parquet")
```

```python
import pyarrow.parquet as pq
study2 = pq.read_table("study2v2.parquet")
```
