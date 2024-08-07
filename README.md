
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PhosphoFind

<!-- badges: start -->
<!-- badges: end -->

The goal of PhosphoFind is to identify the phosphorylation site in a
protein from the phosphosite of a peptide originated from
phosphoproteomics experiments.

## Installation

You can install the development version of PhosphoFind from
[GitHub](https://github.com/RafRomB/PhosphoFind) with:

``` r
# install.packages("devtools")
devtools::install_github("RafRomB/PhosphoFind")
```

## Explanation of PhosphoFind function

The main function in the package is `PhosphoFind`. The input of the
funtion is a dataframe with, at least, the following columns (it is
important that the names of the columns match the ones given here):

- **protein** Uniprot (ACC ID) identifier of the protein.
- **Phospho_sequence**: aminoacid sequence of the peptide.
- **Phosphosite_1**: position in the peptide of the first phosphosite.
- **Phosphosite_2**: position in the peptide of the second phosphosite
  (if any).
- **Pho**: Column indicating if the peptide in that row contains a
  phosphorylation (Y) or not.

The other argument that the function receives, `psp_db`, is a dataframe
with the [PhosphoSitePlus(R)](https://www.phosphosite.org/) [(Hornbeck
et
al. 2014)](https://academic.oup.com/nar/article/43/D1/D512/2439467?login=false)
phosphorylation data for the organism. This dataframe can be loaded with
the function `load_PSP_db`. By default, this function loads the
*Phosphorylation_site_dataset.gz*, Last modified: Fri May 17 09:42:46
EDT 2024, from [PhosphoSitePlus(R)
v6.7.4](https://www.phosphosite.org/staticDownloads) for mouse.
Alternatively, the path to a tab separated value file with a different
PhosphoSitePlus database can be specified through the argument `file`.

### Example

``` r
library(PhosphoFind)

# Load default phosphorylation database
psp_db <- load_PSP_db()

# Choosing a different organism
psp_db <- load_PSP_db(organism = "human")
```

### How it works

The function `PhosphoFind` does the following:

1.  Filters the dataframe to keep only row with a phosphorylated peptide
    `(Pho == "Y")`.
2.  Looks for the protein in the dataframe in the PhosphoSitePlus
    database, based on the Uniprot ACC ID.
3.  If the protein is in the databe, it starts to look for the
    alignments based in the column `SITE_+/-7_AA` from the
    PhosphoSitePlus database (Figure 1) for the first phosphosites in
    the peptides. The column `SITE_+/-7_AA` in the PhosphoSitePlus
    database contains the phosphosite centered by +/- 7 aminoacids (the
    previous 7 aminoacids and the following 7 aminoacids).
4.  If the protein is not in the database, it saves the Uniprot ACC ID
    to return it later.
5.  Looks for aligments in the second phosphosites in the peptides.
6.  Once finished, the function returns a list with two elements. The
    first element is the dataframe with the original dataframe and the
    identified phosphorylation positions in the proteins. The second
    element of the list is a vector with the names of the not identified
    proteins in the database.

<figure>
<img src="man/figures/PhosphoSite.png"
alt="Figure 1. Explanation of ‘PhosphoFind’ function. Depending on the position of the phosphosite in the phosphopeptide, ‘PhosphoFind’ will look for the alignment using different aminoacids. If the phosphosite is at the begining of the phosphopeptide (A), it will use the phosphosite and the following 4 aminoacids. If it is in the second position (B), it will use the previous aminoacid and the following 3 aminoacids to the phosphosite.If it is in the middle (C), it will use the previous 2 and following 2 aminoacids. If it is the second-to-last position (D), it will use the previous 3 and the following aminoacid. And if the phosphosite is at the end of the phosphopepide (E), it will use the previous 4 aminoacids." />
<figcaption aria-hidden="true">Figure 1. Explanation of ‘PhosphoFind’
function. Depending on the position of the phosphosite in the
phosphopeptide, ‘PhosphoFind’ will look for the alignment using
different aminoacids. If the phosphosite is at the begining of the
phosphopeptide (A), it will use the phosphosite and the following 4
aminoacids. If it is in the second position (B), it will use the
previous aminoacid and the following 3 aminoacids to the phosphosite.If
it is in the middle (C), it will use the previous 2 and following 2
aminoacids. If it is the second-to-last position (D), it will use the
previous 3 and the following aminoacid. And if the phosphosite is at the
end of the phosphopepide (E), it will use the previous 4
aminoacids.</figcaption>
</figure>

### Example

``` r
library(PhosphoFind)
data("PhosphoHeart")

# Load 
psp_db <- load_PSP_db(organism = "mouse")

Phospholist <- PhosphoFind(df = PhosphoHeart, psp_db = psp_db)
# Extract dataframe with the identified phosphorylation sites:
Phosphosites <- Phospholist[[1]]
# Extract names of proteins not identified:
No_ID <- Phospholist[[2]]
```

## Additional Functions

### load_PSP_db()

The `load_PSP_db()` function is used to load the phosphorylation sites
information from [PhosphoSitePlus](https://www.phosphosite.org). See
`help(load_PSP_db())` for details.

### phospho_positions()

The `phospho_positions()` function is used to identify the positions of
the phosphosites in a peptide in case these are not available. This
function receives a data frame with the phosphoproteomics data and the
sequence of the modified peptide in a column named ‘peptide’. The
function is specified to detect a particular pattern indicating the
postranslational modifications. Concretely, the function detect the
position of the position of the phosphopeptide if they are specified in
the sequence in the following form:

- *K2(TMT6plex);K33(TMT6plex);S37(Phospho)*

Where `(Phospho)` would indicate the postranslational modification (in
this case, a phosphorylation) and `S37` indicates the aminoacid and the
position in the phosphopeptide of that aminoacid.

In case that the postranslational modifications are indicated in a
different way, another approach should be used to obtain the positions.

# References

Hornbeck PV, Zhang B, Murray B, Kornhauser JM, Latham V, Skrzypek E
[PhosphoSitePlus, 2014: mutations, PTMs and
recalibrations](https://academic.oup.com/nar/article/43/D1/D512/2439467?login=false).
*Nucleic Acids Res*. 2015 43:D512-20. PMID: 25514926.
