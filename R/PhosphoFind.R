#' Find phosphosites
#'
#' This function receives a phosphoproteomics data dataframe with, at least,
#' the following columns:
#' "Pho": Column indicating if the peptides contains a phosphorylation
#' "protein": Uniprot (ACC ID) identifier of the protein
#' "Phospho_sequence": sequence of the peptide
#' "Phosphosite_1": position in the peptide of the first phosphosite
#' "Phosphosite_2": position in the peptide of the second phosphosite (if any)
#'
#' @param df A phosphoproteomics dataframe with the columns previously especified
#' @param psp_db Dataframe with the PhosphoSitePlus phosphorylation data for the
#' organism
#' @return A list with the first element as a dataframe with the phosphoproteomics
#' data with the identified phosphorylation sites. The second element is a vector
#' with the protein Uniprot IDs that have not been identified.
#' @examples
#' data(PhosphoHeart)
#' psp_db <- load_PSP_db(organism = "mouse")
#' # Obtain list with phosphorylation sites and unidentified proteins:
#' Phospholist <- PhosphoFind(df = PhosphoHeart, psp_db = psp_db)
#' # Extract dataframe with the identified phosphorylation sites:
#' Phosphosites <- Phospholist[[1]]
#' # Extract names of proteins not identified:
#' No_ID <- Phospholist[[2]]
#'
#' @export
#'
PhosphoFind <- function(df, psp_db){
  df$ID <- 1:nrow(df)
  df_complete <- df
  df <- df %>% dplyr::filter(Pho == "Y")
  df$Phospho_sequence <-  stringr::str_to_upper(df$Phospho_sequence)
  not_found_IDs <- c()
  pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = nrow(df),
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)
  for ( i in 1:nrow(df)){
    pb$tick()
    if (df$protein[i] %in% psp_db$ACC_ID) { # Check if the protein exists in the database
      phosp_df <- dplyr::inner_join(x = df[i,], y = psp_db, by = dplyr::join_by(protein == ACC_ID)) # Create df with phosphosites
      site_characters <- stringr::str_split(df[i,]$Phospho_sequence, pattern = "", simplify = T) # String to characters
      # If the phosphosite is at the end of the sequence
      if (length(site_characters) == df[i,]$Phosphosite_1) {
        # Characters to string of phosphosite -4 aa
        site_str <- stringr::str_flatten(site_characters[(df[i,]$Phosphosite_1-4):(df[i,]$Phosphosite_1)])
        for (j in 1:nrow(phosp_df)) {
          aa_phosphosite <- stringr::str_split(phosp_df[j,]$"SITE_+/-7_AA", pattern = "", simplify = T)
          aa_str <- stringr::str_flatten(aa_phosphosite[4:8])
          if (aa_str == site_str) {
            df$PHOS_RSD[i] <- phosp_df$MOD_RSD[j]
            break
          } else{
            df$PHOS_RSD[i] <- NA
          }
        }
      }
      # If the phosphosite is at the end -1 of the sequence
      if (length(site_characters) == df[i,]$Phosphosite_1+1) {
        # Characters to string of phosphosite -3 aa and +1 aa
        site_str <- stringr::str_flatten(site_characters[(df[i,]$Phosphosite_1-3):(df[i,]$Phosphosite_1+1)])
        for (j in 1:nrow(phosp_df)) {
          aa_phosphosite <- stringr::str_split(phosp_df[j,]$"SITE_+/-7_AA", pattern = "", simplify = T)
          aa_str <- stringr::str_flatten(aa_phosphosite[5:9])
          if (aa_str == site_str) {
            df$PHOS_RSD[i] <- phosp_df$MOD_RSD[j]
            break
          } else{
            df$PHOS_RSD[i] <- NA
          }
        }
      }
      # If the phosphosite is at the beginning of the sequence
      if (df[i,]$Phosphosite_1 == 1) {
        # Characters to string of phosphosite +4 aa
        site_str <- stringr::str_flatten(site_characters[(df[i,]$Phosphosite_1):(df[i,]$Phosphosite_1+4)])
        for (j in 1:nrow(phosp_df)) {
          aa_phosphosite <- stringr::str_split(phosp_df[j,]$"SITE_+/-7_AA", pattern = "", simplify = T)
          aa_str <- stringr::str_flatten(aa_phosphosite[8:12])
          if (aa_str == site_str) {
            df$PHOS_RSD[i] <- phosp_df$MOD_RSD[j]
            break
          } else{
            df$PHOS_RSD[i] <- NA
          }
        }
      }
      # If the phosphosite is at the +2 position of the sequence
      if (df[i,]$Phosphosite_1 == 2) {
        # Characters to string of phosphosite -1 aa and +3 aa
        site_str <- stringr::str_flatten(site_characters[(df[i,]$Phosphosite_1-1):(df[i,]$Phosphosite_1+3)])
        for (j in 1:nrow(phosp_df)) {
          aa_phosphosite <- stringr::str_split(phosp_df[j,]$"SITE_+/-7_AA", pattern = "", simplify = T)
          aa_str <- stringr::str_flatten(aa_phosphosite[7:11])
          if (aa_str == site_str) {
            df$PHOS_RSD[i] <- phosp_df$MOD_RSD[j]
            break
          } else{
            df$PHOS_RSD[i] <- NA
          }
        }
      }
      # If the phosphosite is in the middle of the sequence
      if(df[i,]$Phosphosite_1 != 1 & length(site_characters) != df[i,]$Phosphosite_1 & length(site_characters) != df[i,]$Phosphosite_1+1 &
         df[i,]$Phosphosite_1 != 2){
        site_str <- stringr::str_flatten(site_characters[(df[i,]$Phosphosite_1-2):(df[i,]$Phosphosite_1+2)])
        for (j in 1:nrow(phosp_df)) {
          aa_phosphosite <- stringr::str_split(phosp_df[j,]$"SITE_+/-7_AA", pattern = "", simplify = T)
          aa_str <- stringr::str_flatten(aa_phosphosite[6:10])
          if (aa_str == site_str) {
            df$PHOS_RSD[i] <- phosp_df$MOD_RSD[j]
            break
          } else{
            df$PHOS_RSD[i] <- NA
          }
        }
      }
    } else{
      #print(paste(df$protein[i],"ID not found"))
      not_found_IDs[i] <- df$protein[i]
      df$PHOS_RSD[i] <- NA
    }
  }
  print("Search for Phosphosites in position 1 finished")
  print("Starting search for position 2")
  pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = nrow(df),
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)
  for (i in 1:nrow(df)){
    pb$tick()
    if(is.na(df$Phosphosite_2[i])){
      df$PHOS_RSD_2[i] <- NA
    }else{
      if (df$protein[i] %in% psp_db$ACC_ID) { # Check if the protein exists in the database
        phosp_df <- dplyr::inner_join(x = df[i,], y = psp_db, by = dplyr::join_by(protein == ACC_ID)) # Create df with phosphosites
        site_characters <- stringr::str_split(df[i,]$Phospho_sequence, pattern = "", simplify = T) # String to characters
        # If the phosphosite is at the end of the sequence
        if (length(site_characters) == df[i,]$Phosphosite_2) {
          # Characters to string of phosphosite -4 aa
          site_str <- stringr::str_flatten(site_characters[(df[i,]$Phosphosite_2-4):(df[i,]$Phosphosite_2)])
          for (j in 1:nrow(phosp_df)) {
            aa_phosphosite <- stringr::str_split(phosp_df[j,]$"SITE_+/-7_AA", pattern = "", simplify = T)
            aa_str <- stringr::str_flatten(aa_phosphosite[4:8])
            if (aa_str == site_str) {
              df$PHOS_RSD_2[i] <- phosp_df$MOD_RSD[j]
              break
            } else{
              df$PHOS_RSD_2[i] <- NA
            }
          }
        }
        # If the phosphosite is at the end -1 of the sequence
        if (length(site_characters) == df[i,]$Phosphosite_2+1) {
          # Characters to string of phosphosite -3 aa and +1 aa
          site_str <- stringr::str_flatten(site_characters[(df[i,]$Phosphosite_2-3):(df[i,]$Phosphosite_2+1)])
          for (j in 1:nrow(phosp_df)) {
            aa_phosphosite <- stringr::str_split(phosp_df[j,]$"SITE_+/-7_AA", pattern = "", simplify = T)
            aa_str <- stringr::str_flatten(aa_phosphosite[5:9])
            if (aa_str == site_str) {
              df$PHOS_RSD_2[i] <- phosp_df$MOD_RSD[j]
              break
            } else{
              df$PHOS_RSD_2[i] <- NA
            }
          }
        }
        # If the phosphosite is at the beginning of the sequence
        if (df[i,]$Phosphosite_2 == 1) {
          # Characters to string of phosphosite +4 aa
          site_str <- stringr::str_flatten(site_characters[(df[i,]$Phosphosite_2):(df[i,]$Phosphosite_2+4)])
          for (j in 1:nrow(phosp_df)) {
            aa_phosphosite <- stringr::str_split(phosp_df[j,]$"SITE_+/-7_AA", pattern = "", simplify = T)
            aa_str <- stringr::str_flatten(aa_phosphosite[8:12])
            if (aa_str == site_str) {
              df$PHOS_RSD_2[i] <- phosp_df$MOD_RSD[j]
              break
            } else{
              df$PHOS_RSD_2[i] <- NA
            }
          }
        }
        # If the phosphosite is at the beginning + 1 of the sequence
        if (df[i,]$Phosphosite_2 == 2) {
          # Characters to string of phosphosite -1 aa and +3 aa
          site_str <- stringr::str_flatten(site_characters[(df[i,]$Phosphosite_2-1):(df[i,]$Phosphosite_2+3)])
          for (j in 1:nrow(phosp_df)) {
            aa_phosphosite <- stringr::str_split(phosp_df[j,]$"SITE_+/-7_AA", pattern = "", simplify = T)
            aa_str <- stringr::str_flatten(aa_phosphosite[7:11])
            if (aa_str == site_str) {
              df$PHOS_RSD_2[i] <- phosp_df$MOD_RSD[j]
              break
            } else{
              df$PHOS_RSD_2[i] <- NA
            }
          }
        }
        # If the phosphosite is in the middle of the sequence
        if(df[i,]$Phosphosite_2 != 1 & df[i,]$Phosphosite_2 != 2 & length(site_characters) != df[i,]$Phosphosite_2 &
           length(site_characters) != df[i,]$Phosphosite_2+1){
          site_str <- stringr::str_flatten(site_characters[(df[i,]$Phosphosite_2-2):(df[i,]$Phosphosite_2+2)])
          for (j in 1:nrow(phosp_df)) {
            aa_phosphosite <- stringr::str_split(phosp_df[j,]$"SITE_+/-7_AA", pattern = "", simplify = T)
            aa_str <- stringr::str_flatten(aa_phosphosite[6:10])
            if (aa_str == site_str) {
              df$PHOS_RSD_2[i] <- phosp_df$MOD_RSD[j]
              break
            } else{
              df$PHOS_RSD_2[i] <- NA
            }
          }
        }
      } else{
        #print(paste(df$protein[i],"ID not found"))
        df$PHOS_RSD_2[i] <- NA
      }
    }
  }
  df <- df %>% dplyr::select(ID, PHOS_RSD, PHOS_RSD_2)
  df <- dplyr::left_join(df_complete, df, by = "ID")
  df <- subset(df, select = -ID)
  return(list("Phosphodf" = df, "IDs" = not_found_IDs))
}
