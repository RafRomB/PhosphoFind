#' Example TMT Phosphoproteomics dataset from Heart Phosphoproteomics
#'
#' Example dataset of a TMT Phosphoproteomics experiment from heart lysates.
#' The dataset corresponds to a edited version of the Supplementary Table 8 from \href{https://www.nature.com/articles/s44161-023-00368-x}{Romero-Becerra et al., 2023}
#' in which the columns with the assigned phosphosites have been omited.
#'
#' Corrected (phosphorylated peptides abundance corrected by protein expression changes) mass spectrometry phosphoproteomics data of hearts
#' lysates from 2-month-old WT mice infected with adeno-associated viruses with cardiac troponin T promoter (AAV-cTnT) either driving expression
#' of constitutively active p38γ/δ kinases (cTnT-p38gamma/delta, n=3) or GFP as a control (cTnT-GFP, n=4). The table includes only phosphopeptides
#' (peptides with at least a phosphorylation as modification).
#' Blank rows indicate phosphopeptides for which only the phosphorylated form was found and normalization was not performed.
#'
#' @docType data
#'
#' @usage data(PhosphoHeart)
#'
#' @format
#' An object of class 'data.frame' with 5327 rows and 16 columns:
#' \describe{
#' \item{**protein**}{Uniprot accession number.}
#' \item{**Uniprot_Entry**}{Uniprot entry name.}
#' \item{**Zpq**}{Indicates phosphopeptide abundance changes corrected by protein expression changes,
#' expressed as log2-ratios in standardized units (Zpq, z-values).}
#' \item{**z_value_difference**}{z values median difference between conditions.}
#' \item{**limma_pvalue**}{p-value associated with the two-sided limma t-test of the z values difference}
#' \item{**Phospho_sequence**}{Sequence of the phosphopeptide.}
#' \item{**Phosphosite_1**}{Position of the first phosphorylated aminoacid in the sequence.}
#' \item{**Phosphosite_2**}{Position of the second phosphorylated aminoacid in the sequence.}
#' }
#'
#' @keywords dataset
#'
#' @source \href{https://www.nature.com/articles/s44161-023-00368-x}{https://www.nature.com/articles/s44161-023-00368-x}
#'
#' @references Romero-Becerra, R., Cruz, F.M., Mora, A. et al. p38γ/δ activation alters cardiac electrical activity and
#' predisposes to ventricular arrhythmia. Nat Cardiovasc Res 2, 1204–1220 (2023). https://doi.org/10.1038/s44161-023-00368-x
#' @examples
#' data(PhosphoHeart)
#' Phospholist <- PhosphoFind(PhosphoHeart)
#' PhosphoHeart <- Phospholist[[1]]
"PhosphoHeart"
