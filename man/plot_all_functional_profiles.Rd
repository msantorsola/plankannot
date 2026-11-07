\name{plot_all_functional_profiles}
\alias{plot_all_functional_profiles}
\title{Generate functional plots for all four databases}
\usage{
plot_all_functional_profiles(
  annotated_files,
  output_dir = "functional_plots",
  trait_map = NULL
)
}
\arguments{
\item{annotated_files}{A named character vector or list of annotation files for the four databases, or a directory produced by \code{\link{annotate_functional_databases}}.}
\item{output_dir}{Directory in which plots are saved.}
\item{trait_map}{Optional trait configuration passed to \code{\link{plot_functional_profiles}}.}
}
\value{
Invisibly returns a named list of results from the database-specific plotting operations.
}
\description{
Resolves the annotation files for the available functional databases and generates database-specific stacked-bar plots and heatmaps using \code{\link{plot_functional_profiles}}.
}
\seealso{
\code{\link{plot_functional_profiles}},
\code{\link{annotate_functional_databases}}
}
