# cloupe

Loupe Browser is a powerful visualization software that provides the intuitive functionality you need to explore and analyze 10x Genomics Chromium and Visium data. This is a PoC parser to read CLOUPE files using python.


# Install

```shell
pip install git+https://github.com/cellgeni/cloupe.git
```

# Run

```python
from cloupe import Cloupe

# parse file
c = Cloupe("/path/to/file.cloupe")

# extract barcodes, features and projections to CSVs
c.to_csv(outdir="/path/to/out")

# create anndata with the same name as input file but `.h5ad` extension
c.to_anndata()
```

### Command line

Export barcodes, features and projections to CSV in current folder
```shell
cloupe --cloupe /path/to/file.cloupe --csv
```