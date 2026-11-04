# REPAIR simulation data

This directory documents the Parquet files used by the REPAIR simulation
studies. The binary datasets are distributed through Zenodo rather than
committed to GitHub: together they are approximately 3.6 GB, and the
individual files exceed GitHub's ordinary repository file-size limit.

## Access

The datasets should be uploaded to a Zenodo data record associated with the
REPAIR software record:

> **REPAIR: Random-effect prediction adjustment for inference in regression**

After Zenodo assigns the data record DOI, replace the placeholder below with
the version-specific DOI and use that DOI when downloading the files:

```text
REPAIR simulation data DOI: <DATA-DOI>
```

The software DOI and data DOI should be cited together when an analysis uses
both the code and these simulation results. The software repository contains
the readers and analysis scripts; this data record contains the generated
simulation tables.

## Files

| File | Size (bytes) | SHA-256 |
|---|---:|---|
| `study1v2.parquet` | 782478809 | `c7c6ef56d6e4265b96a8cf9c007c68b3431ce3cccd00c061cca8384c52463f44` |
| `study2v2.parquet` | 984683117 | `53749e4ffebf73bc6a765220d3551b7b50836f43e2f30da1870982caae485a1f` |
| `study3v2.parquet` | 1920730132 | `ee2db081ba8c40d1e8fb56b7df1803cf4b5f48a3c2532b49c478ef4c12297fab` |
| `study4v2.parquet` | 216399016 | `5aa6e2348f7b40ddf6c4476dcc797233a4142b0d7611a271f33234505d75d162` |

Verify a download with:

```sh
shasum -a 256 study1v2.parquet
shasum -a 256 study2v2.parquet
shasum -a 256 study3v2.parquet
shasum -a 256 study4v2.parquet
```

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
