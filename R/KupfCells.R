#' Example TMT Phosphoproteomics dataset from Kupffer Cells
#'
#' Example dataset of a TMT Phosphoproteomics experiment from Kupffer cells, to
#' exemplify the use of the function 'phospho_positions'. The dataset includes
#' a column with the phosphopeptide sequence and post-translational modifications
#' tags, indluding phosphorylations and others, in the column `peptide`.
#'
#' @docType data
#'
#' @usage data(KupfCells)
#'
#' @format
#' An object of class 'data.frame' with 5000 rows and 6 columns:
#' \describe{
#' \item{**gene**}{Gene symbol.}
#' \item{**protein**}{Uniprot protein ID.}
#' \item{**description**}{Information about the gene encoding the protein,
#' including the gene symbol, complete name, organism, and genome position.}
#' \item{**peptide**}{Column with the peptide and the post-translation modifications.}
#' \item{**Sequence**}{Peptide aminoacid sequence.}
#' \item{**Pho**}{Column indicating if the peptide includes a phosphorylation (`Y`) or not (`NA`)}
#' }
#'
#' @keywords dataset
#'
#' @examples
#' data(KupfCells)
#' KupfCells <- phospho_positions(KupfCells)
"KupfCells"
