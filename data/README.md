# REPAIR simulation data

This directory documents the Parquet files used by the REPAIR simulation
studies. The binary datasets are distributed through Zenodo rather than
committed to GitHub: together they are approximately 3.6 GB, and the
individual files exceed GitHub's ordinary repository file-size limit.

## Access

The datasets should be uploaded to a Zenodo data record associated with the
REPAIR software record:

> **REPAIR: Random-effect prediction adjustment for inference in regression**

The version 1 data record is available at
[10.5281/zenodo.23965547](https://doi.org/10.5281/zenodo.23965547). The
concept DOI for the data across versions is
[10.5281/zenodo.23965546](https://doi.org/10.5281/zenodo.23965546). Use the
version DOI when citing this exact release and the Zenodo record page when
downloading the files:

```text
REPAIR simulation data DOI: 10.5281/zenodo.23965547
REPAIR simulation data concept DOI: 10.5281/zenodo.23965546
```

The software release DOI is
[10.5281/zenodo.23982336](https://doi.org/10.5281/zenodo.23982336). The
software DOI and data DOI should be cited together when an analysis uses both
the code and these simulation results. The software repository contains
the readers and analysis scripts; this data record contains the generated
simulation tables.

## Files

| File | Download | Size (bytes) | MD5 | SHA-256 |
|---|---|---:|---|---|
| `study1v2.parquet` | [Zenodo](https://zenodo.org/records/23965547/files/study1v2.parquet) | 782478809 | `3803f9d7fb432be5dd20f1454c9038c1` | `c7c6ef56d6e4265b96a8cf9c007c68b3431ce3cccd00c061cca8384c52463f44` |
| `study2v2.parquet` | [Zenodo](https://zenodo.org/records/23965547/files/study2v2.parquet) | 984683117 | `f387758cc74cb19c267b885a6b4e01c3` | `53749e4ffebf73bc6a765220d3551b7b50836f43e2f30da1870982caae485a1f` |
| `study3v2.parquet` | [Zenodo](https://zenodo.org/records/23965547/files/study3v2.parquet) | 1920730132 | `8e9b8b05ea846bf7bbab81ff85d82839` | `ee2db081ba8c40d1e8fb56b7df1803cf4b5f48a3c2532b49c478ef4c12297fab` |
| `study4v2.parquet` | [Zenodo](https://zenodo.org/records/23965547/files/study4v2.parquet) | 216399016 | `c778f2c3740ab95061ccb070d7c431bb` | `5aa6e2348f7b40ddf6c4476dcc797233a4142b0d7611a271f33234505d75d162` |

Verify a download with:

```sh
shasum -a 256 study1v2.parquet
shasum -a 256 study2v2.parquet
shasum -a 256 study3v2.parquet
shasum -a 256 study4v2.parquet
```

The MD5 values above are the checksums reported by Zenodo. The SHA-256 values
were computed locally from the source files before upload.

The files are Parquet tables and can be read with R's `arrow` package or
Python's `pyarrow` package. For example:

```r
library(arrow)
study2 <- read_parquet("study2v2.parquet")
```

```python
import pyarrow.parquet as pq
study2 = pq.read_table("study2v2.parquet")
```

The files are release data, not raw participant data. Before publication,
confirm the Zenodo record's metadata and license match the approved data
release.
