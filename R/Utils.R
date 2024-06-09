#' Coalesce across columns.
#'
#' Additional function to coalesce all the columns with the phosphopeptides
#' to remove all the columns with NAs.
#'
#' @export
coacross <- function(...) {
  dplyr::coalesce(!!!dplyr::across(...))
}

#' Find phosphosites positions.
#'
#' This function received a data frame with the phosphoproteomics data and the
#' sequence of the modified peptide in a column named 'peptide' to remove all
#' the undesired modifications that are not a phosphorylation, and identify the
#' position of the phosphorylation site in the phosphopeptide. The function is
#' intended to detect a maximum of two phosphosites per peptide.
#' The input dataframe must have a column call "Pho" wich indicates if the
#' peptide has a phosphorylation in it or not.
#'
#' @param df Dataframe with the phosphoproteomics data.
#' @param ptm Vector with the name of the postranlational modifications in the peptides.
#' different from phosphorylations that are going to be removed.
#' @return A dataframe with the positions of the phosphosites in the peptide.
#' @examples
#' data(KupfCells)
#' phospho_positions(KupfCells)
#'
#' @export
phospho_positions <- function(df, ptm = c("(CARBAMIDOMETHYL)", "(TMT6PLEX)", "(OXIDATION)")){
  df <- as.data.frame(df)
  df$peptide <- stringr::str_to_upper(df$peptide)

  df_ph <- df %>% dplyr::filter(Pho == "Y")

  peptides <- df_ph$peptide %>% stringr::str_split(pattern = ";", simplify = T)

  peptides[,1] <- peptides[,1] %>% stringr::str_replace(pattern = "\\__.*", replacement = "")

  ptm <- ptm

  for (i in 2:ncol(peptides)) {
    for (j in ptm) {
      peptides[,i] <- peptides[,i] %>% stringr::str_replace(pattern = paste(".*",j,".*",sep = ""), replacement = "")
    }
    peptides[,i] <- peptides[,i] %>% stringr::str_replace(pattern = "\\s*\\([^\\)]+\\)", replacement = "")
  }


  df_ph <- cbind(df_ph,peptide_sequence = peptides[,1])

  peptides_df <- as.data.frame(peptides)


  peptides_df[peptides_df == ""] <- NA


  peptides_df <- peptides_df %>% janitor::remove_empty(which = "cols")


  Ph1 <- peptides_df[,2:ncol(peptides_df)] %>% dplyr::mutate(Ph1 = coacross(everything()))

  peptides_df[,2] <- Ph1$Ph1

  for (i in 1:nrow(peptides_df)) {
    for (j in 3:ncol(peptides_df)) {
      if (!is.na(peptides_df[i,j])) {
        if (peptides_df[i,j] == peptides_df[i,2]) {
          peptides_df[i,j] <- NA
        }
      }
    }
  }

  if (ncol(peptides_df)==2) {
    peptides_df$V3 <- NA
  }

  Ph2 <- peptides_df[,3:ncol(peptides_df)] %>% dplyr::mutate(Ph2 = coacross(everything()))
  peptides_df[,3] <- Ph2$Ph2

  peptides_df <- peptides_df[,1:3]


  peptides_df[,2] <- gsub(pattern = "[[:upper:]]", replacement = "", x = peptides_df[,2])
  peptides_df[,3] <- gsub(pattern = "[[:upper:]]", replacement = "", x = peptides_df[,3])


  colnames(peptides_df) <- c("Phospho_sequence", "Phosphosite_1", "Phosphosite_2")

  peptides_df$peptide <- df_ph$peptide


  df <- dplyr::left_join(x = df, y = peptides_df, by = "peptide")

  df$Phosphosite_1 <- as.integer(df$Phosphosite_1)
  df$Phosphosite_2 <- as.integer(df$Phosphosite_2)
  return(df)
}


#' Load phosphosite database.
#'
#' Load phosphorylation sites information from
#' PhosphoSitePlus (https://www.phosphosite.org)
#'
#' @param file Path to the phosphorylation information file. If not especified,
#' the function will use the 'Phosphorylation_site_dataset.gz', Last modified:
#' Fri May 17 09:42:46 EDT 2024 from PhosphoSitePlus(R) v6.7.4 (\href{https://www.phosphosite.org/staticDownloads}{PhosphoSitePlus/Downloads})
#' @param organism String with the name of the organism from which to look the
#' phosphopeptides. Available options so far are: `mouse`, `human`, `rat`, `sheep`,
#' `SARSCoV2`, `guinea pig`, `cow`, `hamster`, `fruit fly`, `dog`, `rabbit`,
#' `pig`, `chicken`, `frog`, `quail`, `horse`, `goat`, `papillomavirus`,
#' `water buffalo`, `marmoset`, `turkey`, `cat`, `starfish`, `torpedo`,
#' `SARSCoV1`, `green monkey`, `ferret`
#' @return A dataframe with the phosphorylation information from the selected.
#' @source The default dataset corresponds to the 'Phosphorylation_site_dataset.gz', Last modified:
#' Fri May 17 09:42:46 EDT 2024 from PhosphoSitePlus(R) v6.7.4 and can be downloaded from
#' \href{https://www.phosphosite.org/staticDownloads}{PhosphoSitePlus/Downloads}
#'
#' PhosphoSite is licensed by Cell Signaling Technology (CST) for non-commercial use under a
#' Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
#' More information about PhosphoSitePlus(R) licensing can be found at \href{https://www.phosphosite.org/staticLicensing.action}{PhosphoSitePlus/About/Licensing}
#' @references Hornbeck PV, Zhang B, Murray B, Kornhauser JM, Latham V, Skrzypek E PhosphoSitePlus, 2014: mutations, PTMs and recalibrations. Nucleic Acids Res. 2015 43:D512-20.
#' @examples
#' # Using the default organims 'mouse'
#' load_PSP_db()
#'
#' # Changing the organism to human
#' load_PSP_db(organism = "human")
#'
#' @export
load_PSP_db <- function(file="PSP_db", organism="mouse"){
  if (file == "PSP_db") {
    phosphosite_db <- PSP_db
  } else{
    phosphosite_db <- readr::read_delim(file = file, delim = "\t")
  }
  phosphosite_db <- phosphosite_db %>% dplyr::filter(ORGANISM == organism)
  phosphosite_db$`SITE_+/-7_AA` <- stringr::str_to_upper(phosphosite_db$`SITE_+/-7_AA`)
  return(phosphosite_db)
}
