#' Parse VIBRANT Quality Score File
#'
#' This function parses a VIBRANT quality score file to extract the quality
#' score and type (lytic or lysogenic) for each contig. It handles common
#' issues in the VIBRANT output, such as duplicate entries and inconsistent
#' classifications. Contig IDs are extracted from scaffold names by removing
#' the "_fragment_*" suffix.
#'
#' @param fin_vibrant_quality A character string specifying the path to
#'   the VIBRANT quality score file.
#'
#' @return A data frame containing the contig ID, quality score, and type
#'   for each contig. Quality scores are re-encoded as numeric values:
#'   * 1: complete circular
#'   * 2: high quality draft
#'   * 3: medium quality draft
#'   * 4: low quality draft
#'   * 9: other (catch-all for unclassified types)
#'
#' @details The function performs the following steps:
#' \enumerate{
#'   \item Reads the VIBRANT quality score file.
#'   \item Extracts the contig ID from scaffold names.
#'   \item Re-encodes the quality score as a numeric value.
#'   \item Separates provirus and non-provirus entries.
#'   \item Resolves duplicate entries, prioritizing non-provirus entries
#'         and higher quality scores.
#'   \item Combines the provirus and non-provirus data.
#'   \item Renames columns to include the "vibrant_" prefix for clarity.
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom data.table fread setnames
#' @importFrom dplyr mutate filter select group_by slice_max slice_min ungroup
#' @importFrom stringr str_replace_all str_detect str_c
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' # Assuming 'vibrant_quality.tsv' is the VIBRANT quality file
#' df_quality <- parse_vibrant_quality("vibrant_quality.tsv")
#' }
parse_vibrant_quality <- function(fin_vibrant_quality) {

  # Step 1: Read the VIBRANT quality file
  df_vibrant <- data.table::fread(fin_vibrant_quality) %>%
    # Extract contig ID from scaffold names by removing "_fragment_*" suffix
    dplyr::mutate(
      contig_id = stringr::str_replace_all(scaffold, "_fragment_.*", "")
    ) %>%
    # Re-encode the quality score as a numeric value using case_when
    dplyr::mutate(
      quality_score = dplyr::case_when(
        Quality == "complete circular" ~ 1,
        Quality == "high quality draft" ~ 2,
        Quality == "medium quality draft" ~ 3,
        Quality == "low quality draft" ~ 4,
        TRUE ~ 9 # Default for unclassified quality levels
      )
    )

  # Step 2: Extract and process provirus entries
  df_vibrant_provirus <- df_vibrant %>%
    # Filter for provirus entries based on scaffold naming pattern
    dplyr::filter(stringr::str_detect(scaffold, "_fragment_")) %>%
    # Select only relevant columns: contig_id and quality_score
    dplyr::select(c("contig_id", "quality_score")) %>%
    # For duplicate contig IDs, keep the entry with the highest quality score
    dplyr::group_by(contig_id) %>%
    dplyr::slice_max(order_by = quality_score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    # Add an empty type column for proviruses
    dplyr::mutate(type = "")

  # Step 3: Extract and process non-provirus entries
  df_vibrant_nonprovirus <- df_vibrant %>%
    # Filter for non-provirus entries
    dplyr::filter(!stringr::str_detect(scaffold, "_fragment_")) %>%
    # Select relevant columns including type
    dplyr::select(c("contig_id", "quality_score", "type")) %>%
    # For duplicate contig IDs, keep the entry with the lowest quality score
    dplyr::group_by(contig_id) %>%
    dplyr::slice_min(order_by = quality_score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  # Step 4: Combine provirus and non-provirus entries
  df_vibrant_quality <- rbind(df_vibrant_nonprovirus, df_vibrant_provirus) %>%
    # Add "vibrant_" prefix to column names for clarity
    data.table::setnames(
      c("quality_score", "type"),
      stringr::str_c("vibrant_", c("quality_score", "type"))
    ) %>%
    # Sort rows by contig_id for easier reference
    dplyr::arrange(contig_id)

  # Return the final processed data frame
  return(df_vibrant_quality)
}


#' Parse VirSorter2 Output
#'
#' This function parses VirSorter2 output files, including the score and
#' category files, to extract relevant information about viral contigs.
#' The function combines information from both files into a single data frame
#' for comprehensive analysis.
#'
#' @param fin_virsorter2_score A character string specifying the path to
#'   the VirSorter2 score file.
#' @param fin_virsorter2_category A character string specifying the path to
#'   the VirSorter2 category file.
#'
#' @return A data frame containing the contig ID, VirSorter2 score, group,
#'   and category for each contig.
#'
#' @details The function performs the following steps:
#' \enumerate{
#'   \item Reads the VirSorter2 score file and renames columns for consistency.
#'   \item Reads the VirSorter2 category file and renames columns.
#'   \item Merges the score and category information on the contig ID.
#'   \item Sorts the final data frame by contig ID for easier reference.
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom data.table fread setnames
#' @importFrom dplyr select inner_join arrange
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' # Assuming 'virsorter2_score.txt' and 'virsorter2_category.txt'
#' # are the VirSorter2 output files
#' df_virsorter2 <- parse_virsorter2(
#'   "virsorter2_score.txt", "virsorter2_category.txt"
#' )
#' }
parse_virsorter2 <- function(fin_virsorter2_score, fin_virsorter2_category) {

  # Step 1: Read the VirSorter2 score file
  df_vs2_summary <- data.table::fread(fin_virsorter2_score) %>%
    # Rename columns for consistency in merged data frame
    data.table::setnames(
      c("seqname", "max_score", "max_score_group"),
      c("contig_id", "virsorter2_score", "virsorter2_group")
    ) %>%
    # Select only relevant columns for analysis
    dplyr::select(
      c("contig_id", "virsorter2_score", "virsorter2_group")
    )

  # Step 2: Read the VirSorter2 category file
  df_vs2_category <- data.table::fread(fin_virsorter2_category) %>%
    # Rename columns for consistent naming conventions
    data.table::setnames(
      c("Contig", "vs2_category"),
      c("contig_id", "virsorter2_category")
    )

  # Step 3: Merge the score and category data frames based on contig ID
  df_vs2 <- df_vs2_summary %>%
    # Inner join ensures only contigs present in both files are included
    dplyr::inner_join(df_vs2_category, by = "contig_id") %>%
    # Sort by contig_id to facilitate ordered data inspection
    dplyr::arrange(contig_id)

  # Return the combined VirSorter2 data frame
  return(df_vs2)
}


#' Parse CheckV Quality
#'
#' This function parses a CheckV quality summary file and extracts
#' relevant information about viral contigs, including their provirus
#' status, quality score, completeness, contamination level, and
#' completeness method. The function re-encodes the `checkv_quality`
#' column as a numeric score for simplified analysis.
#'
#' @param fin_checkv_quality A character string specifying the path to
#'   the CheckV quality summary file.
#'
#' @return A data frame containing the contig ID, provirus status,
#'   re-encoded quality score, completeness, contamination, and the
#'   completeness method for each viral contig.
#'
#' @details The function performs the following steps:
#' \enumerate{
#'   \item Reads the CheckV quality file.
#'   \item Re-encodes the `checkv_quality` column into a numeric quality score.
#'   \item Selects only relevant columns.
#'   \item Renames columns to include the "checkv_" prefix for clarity.
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom data.table fread setnames
#' @importFrom dplyr mutate select case_when
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' # Assuming 'checkv_quality.txt' is the CheckV quality file
#' df_checkv_quality <- parse_checkv_quality("checkv_quality.txt")
#' }
parse_checkv_quality <- function(fin_checkv_quality) {

  # Step 1: Read the CheckV quality file
  df_checkv <- data.table::fread(fin_checkv_quality) %>%
    # Step 2: Re-encode checkv_quality as a numeric quality_score
    dplyr::mutate(
      quality_score = dplyr::case_when(
        checkv_quality == "Complete" ~ 1,
        checkv_quality == "High-quality" ~ 2,
        checkv_quality == "Medium-quality" ~ 3,
        checkv_quality == "Low-quality" ~ 4,
        checkv_quality == "Not-determined" ~ 5,
        TRUE ~ 9 # Default value for unexpected cases
      )
    ) %>%
    # Step 3: Select relevant columns for analysis
    dplyr::select(
      c(
        "contig_id", "contig_length", "provirus", "quality_score",
        "completeness", "contamination", "completeness_method"
      )
    ) %>%
    # Step 4: Rename columns to include a "checkv_" prefix for consistency
    data.table::setnames(
      c(
        "contig_id", "contig_length", "checkv_provirus", "checkv_quality_score",
        "checkv_completeness", "checkv_contamination",
        "checkv_completeness_method"
      )
    ) %>%
    dplyr::arrange(contig_id)

  # Return the parsed CheckV quality data frame
  return(df_checkv)
}



#' Parse AniCLST Output (Clusters based on Average Nucleic Acid Identity)
#'
#' This function parses the AniCLST output, which maps contigs to their
#' corresponding representative contig IDs (Rep IDs). The function assigns
#' unique vOTU IDs to each Rep ID and creates a data frame mapping each
#' contig ID to a vOTU ID.
#'
#' @param fin_aniclst A character string specifying the path to the AniCLST
#'   output file.
#'
#' @return A data frame with columns `contig_repid`, `contig_id` and `vOTU_id`,
#'   which maps contig IDs to representative contig IDs and unique vOTU IDs for
#'   easy identification of clusters.
#'
#' @details The function performs the following steps:
#' \enumerate{
#'   \item Reads the AniCLST output file.
#'   \item Renames columns for consistent naming conventions.
#'   \item Generates a unique `vOTU_id` for each unique Rep ID.
#'   \item Maps `contig_id` values to `vOTU_id` values for downstream use.
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom data.table fread setnames
#' @importFrom dplyr arrange select distinct mutate left_join
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' # Assuming 'aniclst_output.txt' is the AniCLST output file
#' df_ctg2vOTU <- parse_aniclst("aniclst_output.txt")
#' }
parse_aniclst <- function(fin_aniclst) {
  # Step 1: Read the AniCLST output file
  df_aniclst <- data.table::fread(fin_aniclst) %>%
    # Rename columns for consistency with internal naming conventions
    data.table::setnames(c("contig_repid", "contig_id"))

  # Step 2: Generate unique vOTU IDs for each Rep ID
  df_repid <- df_aniclst %>%
    # Arrange by Rep ID for consistent vOTU ID assignment order
    dplyr::arrange(contig_repid) %>%
    # Select and keep only unique Rep IDs
    dplyr::select("contig_repid") %>%
    dplyr::distinct() %>%
    # Create unique vOTU IDs with a zero-padded format (e.g., vOTU000001)
    dplyr::mutate(
      n_digits = ceiling(log10(n())),
      vOTU_id = paste0("vOTU", sprintf(paste0("%0", n_digits, "d"), 1:n()))
    ) %>%
    dplyr::select(-c("n_digits"))

  # Step 3: Map contig IDs to their corresponding vOTU IDs
  df_ctg2vOTU <- df_aniclst %>%
    # Perform a left join to assign each contig_id the correct vOTU_id
    dplyr::left_join(df_repid, by = "contig_repid") %>%
    # Arrange by vOTU_id for organized output
    dplyr::arrange(vOTU_id)

  # Return the final data frame mapping contig IDs to vOTU IDs
  return(df_ctg2vOTU)
}


# Stable ties preserve the member order used by the submitted catalog. Animal
# source fields summarize the cluster; other fields describe the selected row.
build_votu_representatives <- function(filtered, members) {
  sources <- members %>%
    dplyr::select(vOTU_id, dplyr::starts_with("host_")) %>%
    dplyr::group_by(vOTU_id) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop")
  filtered %>%
    dplyr::arrange(dplyr::desc(checkv_completeness), contig_id) %>%
    dplyr::distinct(vOTU_id, .keep_all = TRUE) %>%
    dplyr::mutate(original_representative_id = contig_repid,
      contig_repid = contig_id) %>%
    dplyr::select(-dplyr::starts_with("host_")) %>%
    dplyr::left_join(sources, by = "vOTU_id", relationship = "one-to-one") %>%
    dplyr::mutate(vOTU_prevalence = rowSums(dplyr::pick(dplyr::starts_with("host_")) > 0)) %>%
    dplyr::select(-dplyr::all_of(c("contig_id", "study_id", "assembly_id", "sample_id")))
}


#' Parse vContact3 Output
#'
#' This function parses the vContact3 output file, specifically the
#' `final_assignments.csv` file. It processes genus predictions,
#' classifies contigs, and identifies novel genera with or without
#' reference genomes.
#'
#' @param fin_vcontact3 Character string specifying the path to the
#'   vContact3 `final_assignments.csv` file.
#'
#' @return A data frame containing parsed vContact3 information with
#'   VC IDs, classification status, reference presence, and novel
#'   genus flags.
#'
#' @importFrom magrittr %>%
#' @importFrom data.table fread
#' @importFrom data.table setnames
#' @importFrom dplyr mutate
#' @importFrom dplyr if_else
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr distinct
#' @importFrom dplyr left_join
#' @importFrom stringr str_detect
#'
#' @examples
#' \dontrun{
#' # Assuming 'final_assignments.csv' is the vContact3 output file
#' df_vcontact3 <- parse_vcontact3("final_assignments.csv")
#' }
parse_vcontact3 <- function(fin_vcontact3) {
  # Read the vContact3 final_assignments.csv file.
  df_vc3 <- data.table::fread(fin_vcontact3) %>%
    # Rename the prediction column for clarity.
    data.table::setnames(c("genus (prediction)"), c("genus_prediction")) %>%
    # If genus_prediction is empty, use GenomeName as a placeholder.
    # This is because GenomeName is usually a species (when using vOTU
    # representative contigs for vConTACT3), not a genus, but we have no
    # other choice. Using GenomeName is reasonable because these contigs
    # do not resemble other contigs or viruses in both vConTACT3 and
    # reference genomes. This indicates they likely belong to a novel
    # genus, or at least a novel species. Additionally, such contigs
    # often have short genome lengths (<10 kb), indicating low completeness.
    dplyr::mutate(
      vc_genus = dplyr::if_else(
        genus_prediction == "", GenomeName, genus_prediction
      )
    ) %>%
    # Select relevant columns.
    dplyr::select("GenomeName", "vc_genus", "Reference") %>%
    # Classify contigs based on whether a genus was assigned. Here,
    # "classified" means that vConTACT3 assigns a genus to the contig.
    # The genus name could be either a true virus genus or a novel genus ID.
    # "unclassified" indicates that vConTACT3 cannot assign a genus to the
    # contig, so we use GenomeName instead as a placeholder.
    dplyr::mutate(
      vc_classified = dplyr::if_else(
        vc_genus == GenomeName, "unclassified", "classified"
      )
    )

  # Generate unique VC IDs for each vc_genus.
  # Warning: The VC IDs created here include reference genomes, so the total
  # number of VCs does not represent only query genomes. We'll calculate
  # VCs containing both reference and query genomes later if needed.
  vc_ids <- df_vc3 %>%
    dplyr::select("vc_genus") %>%
    dplyr::distinct() %>%
    # Create unique VC IDs with a zero-padded format (e.g., VC000001)
    # for each unique genus.
    dplyr::mutate(
      n_digits = ceiling(log10(n())),
      vc_id = paste0("VC", sprintf(paste0("%0", n_digits, "d"), 1:n()))
    ) %>%
    # Drop the n_digits column, as it is no longer needed.
    dplyr::select(-"n_digits")

  # Join the VC IDs with the main data frame.
  df_vc3_with_vcid <- df_vc3 %>%
    dplyr::left_join(vc_ids, by = "vc_genus")

  # Identify VCs that contain reference genomes.
  vcs_with_ref <- df_vc3_with_vcid %>%
    dplyr::filter(Reference == "TRUE") %>%
    dplyr::select("vc_id") %>%
    dplyr::distinct()

  # Add a flag for VCs with references and identify novel genera.
  df_vc3_with_vcid <- df_vc3_with_vcid %>%
    dplyr::mutate(
      # Flag VCs that include reference genomes.
      vc_with_ref = dplyr::if_else(
        vc_id %in% vcs_with_ref$vc_id, "TRUE", "FALSE"
      ),
      # Flag as novel genera if a VC doesn't contain reference genomes
      # and the genus name starts with "novel_". Note that this excludes
      # "unclassified" query contigs since we lack genus predictions for them.
      vc_novel_genus = dplyr::if_else(
        stringr::str_detect(vc_genus, "^novel_") & vc_with_ref == "FALSE",
        "TRUE", "FALSE"
      )
    )

  # Return the parsed vContact3 data frame with VC IDs and classification.
  return(df_vc3_with_vcid)
}



#' Parse BACPHLIP Output
#'
#' This function parses the output of BACPHLIP software, which
#' predicts the lifestyle (temperate or virulent) of phage genomes.
#'
#' @param fin_bacphlip Character string specifying the path to the
#'   BACPHLIP output file.
#'
#' @return A data frame containing contig IDs and predicted
#'   lifestyles as `bacphlip_lifestyle`.
#'
#' @importFrom magrittr %>%
#' @importFrom data.table fread
#' @importFrom dplyr select
#' @importFrom data.table setnames
#'
#' @examples
#' \dontrun{
#' # Assuming 'bacphlip_output.txt' is the BACPHLIP output file
#' df_bacphlip <- parse_bacphlip("bacphlip_output.txt")
#' }
parse_bacphlip <- function(fin_bacphlip) {
  # Read the BACPHLIP output file.
  df <- data.table::fread(fin_bacphlip) %>%
    # Select relevant columns: "Contig" (phage genome ID) and "replicyc"
    # (predicted lifestyle).
    dplyr::select("Contig", "replicyc") %>%
    # Rename columns for clarity: "Contig" to "contig_id" and
    # "replicyc" to "bacphlip_lifestyle".
    data.table::setnames(c("contig_id", "bacphlip_lifestyle"))

  # Return the parsed data frame with contig IDs and predicted lifestyle.
  return(df)
}


parse_genomad_summary <- function(path_genomad_summary) {
  df <- data.table::fread(path_genomad_summary) %>%
    # split the taxonomy string into 7 columns
    data.table::setnames(colnames(.), str_c("genomad_", colnames(.))) %>%
    separate(genomad_taxonomy, into = c("classify_type", "realm", "kingdom", "phylum", "class", "order", "family"), sep = ";", fill = "right") %>%
    dplyr::mutate(contig_id = str_replace_all(genomad_seq_name, "\\|.*", "")) %>%
    dplyr::select(-c("genomad_seq_name", "genomad_length", "genomad_coordinates")) %>%
    # put `ctgid` as the first column
    dplyr::select(contig_id, everything()) %>%
    dplyr::mutate(genus = "", species = "")

  return(df)
}


parse_genomad_taxonomy <- function(df_genomad) {
  df_genomad_taxonomy <-df_genomad %>%
    dplyr::filter(classify_type == "Viruses") %>%
    dplyr::select(contig_id, realm, kingdom, phylum, class, order, family, genus, species) %>%
    # remove rows where all taxonomy columns are empty
    dplyr::filter(!dplyr::if_all(c("realm", "kingdom", "phylum", "class", "order", "family", "genus", "species"), ~ . == "")) %>%
    # fill taxonomy columns with "[Taxa]_<Taxa>" if the column is empty but the next column is not empty
    dplyr::mutate(
      genus   = if_else(genus == "", paste0("[Genus]_", species), genus),
      family  = if_else(family == "", paste0("[Family]_", str_replace_all(genus, "\\[Genus\\]_", "")), family),
      order   = if_else(order == "", paste0("[Order]_", str_replace_all(family, "\\[Family\\]_", "")), order),
      class   = if_else(class == "", paste0("[Class]_", str_replace_all(order, "\\[Order\\]_", "")), class),
      phylum  = if_else(phylum == "", paste0("[Phylum]_", str_replace_all(class, "\\[Class\\]_", "")), phylum),
      kingdom = if_else(kingdom == "", paste0("[Kingdom]_", str_replace_all(phylum, "\\[Phylum\\]_", "")), kingdom),
      realm   = if_else(realm == "", paste0("[Realm]_", str_replace_all(kingdom, "\\[Kingdom\\]_", "")), realm),
      # # remove cells ending with "_" (i.e. incomplete taxonomic names)
      across(everything(), ~ if_else(str_ends(., "_"), "", .))
    )

  return(df_genomad_taxonomy)
}


parse_vitap <- function(path_vitap, contig_marker = "ctg") {
  df_vitap <- data.table::fread(path_vitap) %>%
    dplyr::filter(str_detect(Genome_ID, contig_marker)) %>%
    tidyr::separate(lineage, into = c("species", "genus", "family", "order", "class", "phylum", "kingdom", "realm"), sep = ";") %>%
    dplyr::select(Genome_ID, realm, kingdom, phylum, class, order, family, genus, species, lineage_score, Confidence_level) %>%
    data.table::setnames(c("Genome_ID", "Confidence_level"), c("contig_id", "lineage_confidence")) %>%
    dplyr::select(-c("lineage_score", "lineage_confidence")) %>%
    # replace all cells value with "-" to ""
    dplyr::mutate(across(everything(), ~ ifelse(. == "-", "", .)))

  return(df_vitap)
}


read_amg <- function(fin_amg, amg_source) {
  if (amg_source == "vibrant") {
    df_amg <- data.table::fread(fin_amg) %>%
      data.table::setnames(c("protein", "AMG KO", "AMG KO name"), c("protein_id", "dbid", "db_desc")) %>%
      dplyr::mutate(contig_id = str_replace_all(scaffold, "_fragment_[0-9]+$", ""))
  } else if (amg_source == "dramv") {
    df_amg <- data.table::fread(fin_amg) %>%
      dplyr::mutate(contig_id = str_replace_all(scaffold, "-cat_[0-9]+$", "")) %>%
      data.table::setnames(c("gene_id", "gene_description", "gene"), c("dbid", "db_desc", "orf_id")) %>%
      dplyr::mutate(protein_id = str_replace_all(orf_id, "-cat_[0-9]+_", "_"))
  } else {
    stop("Invalid AMG source")
  }

  df_amg <- df_amg %>%
    dplyr::select(contig_id, protein_id, dbid, db_desc, everything()) %>%
    data.table::setnames(colnames(.), make.names(colnames(.)))

  return(df_amg)
}




#' Convert a taxonomy data frame to a Sankey plot data frame
#'
#' This function converts a taxonomy data frame into a Sankey plot data frame,
#' which can be used to create a Sankey plot using the `networkD3` package.
#'
#' @param df_taxonomy A data frame containing taxonomy information with
#'   columns representing different taxonomic levels.
#'
#' @param level_cols A character vector specifying the taxonomic levels
#'   to include in the Sankey plot.
#'
#' @return A list containing two data frames:
#'   - `links`: A data frame with source and target columns for the Sankey plot.
#'   - `nodes`: A data frame with node names and group assignments for coloring nodes.
#'
#' @importFrom magrittr %>%
#' @importFrom data.table fread
#' @importFrom data.table setnames
#' @importFrom dplyr mutate
#' @importFrom dplyr if_else
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr distinct
#' @importFrom dplyr left_join
#' @importFrom stringr str_detect
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr rename
#' @importFrom dplyr sym
#' @importFrom dplyr all_of
#' @importFrom networkD3 sankeyNetwork
#' @importFrom tidyr separate
#' @importFrom tidyr unite
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr fill
#' @importFrom tidyr replace_na
#'
#' @examples
#' \dontrun{
#' # Assuming 'df_taxonomy' is the taxonomy data frame
#' df_sankey <- df_to_sankey(df_taxonomy, c("realm", "kingdom", "phylum", "class", "order", "family", "genus", "species"))
#' }
df_to_sankey <- function(df_taxonomy, level_cols) {
  # --- 1. Prepare Links Data ---
  # This section transforms your table into a list of source-target links with counts.
  # Initialize an empty list to store the links
  links_list <- list()
  for (i in 1:(length(level_cols) - 1)) {
    # Check if both source and target columns exist in the dataframe
    source_col <- level_cols[i]
    target_col <- level_cols[i+1]

    if (!(source_col %in% names(df_taxonomy)) || !(target_col %in% names(df_taxonomy))) {
      warning(paste("Skipping link creation between", source_col, "and", target_col, "- one or both columns not found."))
      next # Skip to the next iteration if columns are missing
    }

    # --- MODIFIED SECTION START ---
    link_data <- df_taxonomy %>%
      # Select the source and target columns using their string names safely
      dplyr::select(all_of(c(source_col, target_col))) %>%
      # Ensure columns are character type using across()
      dplyr::mutate(across(all_of(c(source_col, target_col)), as.character)) %>%
      # Remove rows where source or target might be NA or empty strings
      # Use .data[[var]] notation for filtering with string variables
      dplyr::filter(!is.na(.data[[source_col]]), !is.na(.data[[target_col]]),
             .data[[source_col]] != "", .data[[target_col]] != "") %>%
      # Count the occurrences of each unique link (source-target pair)
      # Group by the original column names using !!sym()
      dplyr::group_by(!!dplyr::sym(source_col), !!dplyr::sym(target_col)) %>%
      dplyr::summarise(value = dplyr::n(), .groups = 'drop') %>%
      # Rename the columns to 'source' and 'target' AFTER summarising
      dplyr::rename(source = !!dplyr::sym(source_col), target = !!dplyr::sym(target_col))
    # --- MODIFIED SECTION END ---

    if (nrow(link_data) > 0) {
        links_list[[length(links_list) + 1]] <- link_data
    }
  }

  # Combine all links from different levels into one data frame
  if (length(links_list) > 0) {
    links <- dplyr::bind_rows(links_list)
  } else {
    stop("Error: No valid links could be created. Check column names and data content.")
  }


  # --- 2. Create Nodes Data Frame ---

  # Get unique node names from all valid sources and targets
  nodes <- data.frame(name = unique(c(links$source, links$target)), stringsAsFactors = FALSE)

  # Add a 'group' column (optional) - for coloring nodes by taxonomic rank
  # This assigns the group based on the highest taxonomic rank the name appears in
  nodes$group <- NA # Initialize
  for(col_name in level_cols) {
    if (col_name %in% names(df_taxonomy)) { # Check if column exists before accessing
        # Ensure column exists before trying to access it
        if(col_name %in% names(df_taxonomy)) {
            current_level_nodes <- unique(df_taxonomy[[col_name]])
            # Ensure comparison is done with character type
            current_level_nodes <- as.character(current_level_nodes)
            nodes$group[nodes$name %in% current_level_nodes & is.na(nodes$group)] <- col_name
        }
    }
  }
  # Assign a default group if any NA remains
  nodes$group[is.na(nodes$group)] <- "Unknown"


  # --- 3. Map Node Names to Zero-Based IDs ---
  # networkD3 requires numeric, 0-indexed IDs for sources and targets.

  # Create the 0-based ID for each node
  nodes$ID <- 0:(nrow(nodes) - 1)

  # Add source and target IDs to the links data frame by matching names
  # Ensure links$source and links$target are character for matching
  links$source <- as.character(links$source)
  links$target <- as.character(links$target)
  nodes$name <- as.character(nodes$name)

  links$sourceID <- match(links$source, nodes$name) - 1
  links$targetID <- match(links$target, nodes$name) - 1

  # Check for NAs in IDs - indicates a mismatch if any exist
  if(any(is.na(links$sourceID)) || any(is.na(links$targetID))) {
      warning("NA values generated during ID matching. Check consistency between links and nodes.")
      # Optional: Filter out links with NA IDs if needed
      # links <- links %>% filter(!is.na(sourceID), !is.na(targetID))
  }


  # Prepare the final links data frame for plotting (only needs IDs and value)
  # Ensure source/target are integers
  links_for_plot <- links %>%
      dplyr::filter(!is.na(sourceID), !is.na(targetID)) %>% # Remove links with potential NA IDs
      dplyr::select(source = sourceID, target = targetID, value) %>%
      dplyr::mutate(source = as.integer(source), target = as.integer(target)) %>%
      # convert to a plain data frame
      as.data.frame()

  return(list(links = links_for_plot, nodes = nodes))
}
