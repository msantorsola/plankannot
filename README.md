# plankannot

`plankannot` annotates user-provided plankton taxonomy tables against four
curated functional databases (Copepoda, HABs, Mixoplankton and Phytoplankton)
and generates abundance-based functional visualizations.

## Installation

The package is intended for distribution through CRAN:

```r
install.packages("plankannot")
library(plankannot)
```

During development, install the source package from a local checkout with:

```r
# install.packages("devtools")
devtools::install()
```

## Input data

`annotate_functional_databases()` accepts XLSX, CSV, TSV and TXT files. The
input must contain at least two columns. The first column is retained as the
feature or sequence identifier. Values in subsequent columns are scanned for
candidate taxonomy or species names. Semicolon-delimited taxonomy strings are
split into individual taxon names.

The function compares candidate names with the `Species_Name` field in each
selected bundled database. Matching is based on normalized text values; it is
not fuzzy taxonomic identification. The function does not infer taxonomic
assignments.

For plotting, sample-abundance columns should be named `Sample_1`,
`Sample_2`, and so on. The plotting functions detect these columns by name.

## Workflow

### 1. Annotate against the functional databases

```r
library(plankannot)

annotation_results <- annotate_functional_databases(
  input_file = "path/to/your_input.xlsx",
  output_dir = "functional_analysis/annotations"
)
```

The function writes database-specific TSV files and, by default, XLSX files,
as well as a CSV matching summary. The output tables retain the input columns
and include matching information.

To process only selected databases, use `db_names`:

```r
annotate_functional_databases(
  input_file = "path/to/your_input.xlsx",
  db_names = c("copepoda", "habs"),
  output_dir = "functional_analysis/annotations"
)
```

### 2. Generate functional plots

```r
plot_all_functional_profiles(
  annotated_files = "functional_analysis/annotations",
  output_dir = "functional_analysis/plots"
)
```

Alternatively, generate plots for one database:

```r
plot_functional_profiles(
  annotated_data = "functional_analysis/annotations/path_to_annotation.xlsx",
  db_name = "copepoda",
  output_dir = "functional_analysis/plots"
)
```

The package includes predefined functional traits for each database. The
available traits are determined by the bundled database annotations.

## Bundled functional databases

The package includes four RDS databases in `inst/db/`:

- `copepoda.rds`
- `habs.rds`
- `mixoplankton.rds`
- `phytoplankton.rds`

See `inst/DATA_SOURCES.md` for information about their preparation and
provenance.

## Citation

Use `citation("plankannot")` to obtain the package citation.

## Authors and contributions

The functional databases and package were developed by Said Ahmed and
Mariangela Santorsola.

## Reporting issues

Please report bugs and feature requests at the package issue tracker:
https://github.com/msantorsola/plankannot/issues
