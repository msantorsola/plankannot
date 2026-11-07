\name{plot_functional_profiles}
\alias{plot_functional_profiles}
\title{Generate functional stacked-bar plots and heatmaps}
\usage{
plot_functional_profiles(
  annotated_data,
  db_name,
  output_dir = "functional_plots",
  sample_cols = NULL,
  trait_map = NULL
)
}
\arguments{
\item{annotated_data}{A data frame or path to an XLSX or TSV annotation file.}
\item{db_name}{Name of the database represented by the annotation data. One of \code{"mixoplankton"}, \code{"copepoda"}, \code{"phytoplankton"}, or \code{"habs"}.}
\item{output_dir}{Directory in which plots are saved.}
\item{sample_cols}{Optional character vector of sample-abundance column names. If \code{NULL}, columns named \code{Sample_N} are detected.}
\item{trait_map}{Optional named character vector or list defining the functional traits to plot.}
}
\value{
A list containing the detected traits, plot data, and output directory information.
}
\description{
Reads an annotation table and summarizes selected categorical functional traits using sample abundances. For each available selected trait, creates a stacked-bar plot and a heatmap.
}
\seealso{
\code{\link{annotate_functional_databases}},
\code{\link{plot_all_functional_profiles}}
}
