% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utils.R
\name{load_PSP_db}
\alias{load_PSP_db}
\title{Load phosphosite database.}
\source{
The default dataset corresponds to the 'Phosphorylation_site_dataset.gz', Last modified:
Fri May 17 09:42:46 EDT 2024 from PhosphoSitePlus(R) v6.7.4 and can be downloaded from
\href{https://www.phosphosite.org/staticDownloads}{PhosphoSitePlus/Downloads}

PhosphoSite is licensed by Cell Signaling Technology (CST) for non-commercial use under a
Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
More information about PhosphoSitePlus(R) licensing can be found at \href{https://www.phosphosite.org/staticLicensing.action}{PhosphoSitePlus/About/Licensing}
}
\usage{
load_PSP_db(file = "PSP_db", organism = "mouse")
}
\arguments{
\item{file}{Path to the phosphorylation information file. If not especified,
the function will use the 'Phosphorylation_site_dataset.gz', Last modified:
Fri May 17 09:42:46 EDT 2024 from PhosphoSitePlus(R) v6.7.4 (\href{https://www.phosphosite.org/staticDownloads}{PhosphoSitePlus/Downloads})}

\item{organism}{String with the name of the organism from which to look the
phosphopeptides. Available options so far are: \code{mouse}, \code{human}, \code{rat}, \code{sheep},
\code{SARSCoV2}, \verb{guinea pig}, \code{cow}, \code{hamster}, \verb{fruit fly}, \code{dog}, \code{rabbit},
\code{pig}, \code{chicken}, \code{frog}, \code{quail}, \code{horse}, \code{goat}, \code{papillomavirus},
\verb{water buffalo}, \code{marmoset}, \code{turkey}, \code{cat}, \code{starfish}, \code{torpedo},
\code{SARSCoV1}, \verb{green monkey}, \code{ferret}}
}
\value{
A dataframe with the phosphorylation information from the selected.
}
\description{
Load phosphorylation sites information from
PhosphoSitePlus (https://www.phosphosite.org)
}
\examples{
# Using the default organims 'mouse'
load_PSP_db()

# Changing the organism to human
load_PSP_db(organism = "human")

}
\references{
Hornbeck PV, Zhang B, Murray B, Kornhauser JM, Latham V, Skrzypek E PhosphoSitePlus, 2014: mutations, PTMs and recalibrations. Nucleic Acids Res. 2015 43:D512-20.
}
