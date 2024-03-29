% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peptideLevel_DE.R
\name{plot_1prot}
\alias{plot_1prot}
\title{Plot trends for a single protien}
\usage{
plot_1prot(
  mm,
  prot.info,
  prot_to_plot,
  prot_to_plot_col,
  gene_name,
  gene_name_col,
  colors,
  mylabs
)
}
\arguments{
\item{mm}{matrix of raw intensities}

\item{prot.info}{metadata for the intensities in mm}

\item{prot_to_plot}{protein ID to plot}

\item{prot_to_plot_col}{protein ID column index}

\item{gene_name}{gene ID to plot}

\item{gene_name_col}{gene ID to plot column index}

\item{colors}{what colors to plot peptide abundances as,
most commonly should be
treatment groups}

\item{mylabs}{sample labels to be plotted on x-axis}
}
\value{
Nil
}
\description{
Plot peptide trends for a protein
}
