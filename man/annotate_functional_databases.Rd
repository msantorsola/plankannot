\name{annotate_functional_databases}
\alias{annotate_functional_databases}
\title{Annotate an input table against the bundled functional databases}
\usage{
annotate_functional_databases(
  input_file,
  output_dir = "functional_analysis/annotations",
  db_names = c("copepoda", "mixoplankton", "phytoplankton", "habs"),
  db_dir = NULL,
  write_xlsx = TRUE
)
}
\arguments{
\item{input_file}{Path to an XLSX, CSV, TSV, or TXT input file. The input must contain at least two columns. The first column is retained as the feature or sequence identifier; values in subsequent columns are scanned for taxonomy or species names. Semicolon-delimited taxonomy strings are split into individual candidate names.}
\item{output_dir}{Directory in which database-specific annotation results are written.}
\item{db_names}{Character vector specifying which bundled databases to use. Available names are \code{"copepoda"}, \code{"mixoplankton"}, \code{"phytoplankton"}, and \code{"habs"}.}
\item{db_dir}{Directory containing the database RDS files. If \code{NULL}, the package database directory is used.}
\item{write_xlsx}{Logical; whether to write XLSX files in addition to TSV files.}
}
\value{
Invisibly returns a list containing output file paths, matching summaries, detected sample columns, and scanned taxonomy columns.
}
\description{
Scans input values after the first column and matches candidate taxon names against the \code{Species_Name} field in the selected functional databases. Matching is based on normalized text values; semicolon-delimited taxonomy strings are split and scanned from right to left. Writes database-specific annotation tables and a matching summary. The function does not perform taxonomic identification.
}
\details{
The input file may be in XLSX, CSV, TSV, or TXT format. For CSV files, comma is used as the field delimiter; for TSV and TXT files, tab is used. Matching is based on names present in the input values and the database \code{Species_Name} column. The function does not perform fuzzy matching.
}
\seealso{
\code{\link{plot_functional_profiles}},
\code{\link{plot_all_functional_profiles}}
}
