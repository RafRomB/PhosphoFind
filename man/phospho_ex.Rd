% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phospho_ex.R
\docType{data}
\name{phospho_ex}
\alias{phospho_ex}
\title{Example TMT Phosphoproteomics dataset}
\format{
An object of class 'data.frame' with 5000 rows and 6 columns:
\describe{
\item{\strong{gene}}{Gene symbol.}
\item{\strong{protein}}{Uniprot protein ID.}
\item{\strong{description}}{Information about the gene encoding the protein,
including the gene symbol, complete name, organism, and genome position.}
\item{\strong{peptide}}{Column with the peptide and the post-translation modifications.}
\item{\strong{Sequence}}{Peptide aminoacid sequence.}
\item{\strong{Pho}}{Column indicating if the peptide includes a phosphorylation (\code{Y}) or not (\code{NA})}
}
}
\usage{
data(phospho_ex)
}
\description{
Example dataset of a TMT Phosphoproteomics experiment, to
exemplify the use of the function 'phospho_positions'. The dataset includes
a column with the phosphopeptide sequence and post-translational modifications
tags, indluding phosphorylations and others, in the column \code{peptide}.
}
\examples{
data(phospho_ex)
phospho_ex <- phospho_positions(phospho_ex)
}
\keyword{dataset}
