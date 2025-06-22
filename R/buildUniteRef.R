# Biostrings required for fasta file handling
library(Biostrings)

#' Create formatted reference fasta and taxonomy files from a UNITE general release
#'
#' The buildUniteRef function processes a UNITE reference fasta file by extracting and
#' filtering taxonomic annotations. It validates species names, handles problematic labels
#' and produces two separate output files: a reference fasta sequence file and a
#' tab-separated file containing reference ID and filtered taxonomic classification.
#'
#'
#' @importFrom Biostrings DNAStringSet readDNAStringSet readRNAStringSet writeXStringSet
#' @importFrom utils read.table write.table
#'
#' @param fin Path to the input UNITE FASTA file.
#' @param fout_fasta Path for the output FASTA file.
#' @param fout_taxonomy Path for the output taxonomy TSV file.
#' @param include.species Logical; if TRUE, includes validated species names.
#' @param compress Logical; if TRUE, the output FASTA file will be gzip compressed.
#'
#' @return Function writes two files to disk, but does not return values.
#' @export
#'
#' @examples
#' \dontrun{
#' buildUniteRef("sh_general_release_dynamic_10.05.2021.fasta",
#'               "unite_ref.fasta", "unite_taxonomy.tsv",
#'               include.species = TRUE, compress = FALSE)
#' }
buildUniteRef <- function(fin, fout_fasta, fout_taxonomy,
                          include.species = TRUE,
                          compress = FALSE) {

    # Start timing
    start_time <- Sys.time()

    # Parameter validation
    if (missing(fin) || !is.character(fin) || length(fin) != 1) {
        stop("'fin' must be a single character string specifying the input FASTA file path.")
    }
    if (missing(fout_fasta) || !is.character(fout_fasta) || length(fout_fasta) != 1) {
        stop("'fout_fasta' must be a single character string specifying the output FASTA file path.")
    }
    if (missing(fout_taxonomy) || !is.character(fout_taxonomy) || length(fout_taxonomy) != 1) {
        stop("'fout_taxonomy' must be a single character string specifying the output taxonomy file path.")
    }
    if (!is.logical(include.species) || length(include.species) != 1) {
        stop("'include.species' must be a single logical value (TRUE or FALSE).")
    }
    if (!is.logical(compress) || length(compress) != 1) {
        stop("'compress' must be a single logical value (TRUE or FALSE).")
    }
    if (!file.exists(fin)) {
        stop("Input FASTA file not found: ", fin)
    }

    # Create output directories if they don't exist
    fout_fasta_dir <- dirname(fout_fasta)
    fout_taxonomy_dir <- dirname(fout_taxonomy)
    if (!dir.exists(fout_fasta_dir)) {
        tryCatch({
            dir.create(fout_fasta_dir, recursive = TRUE, showWarnings = FALSE)
            cat("Created output directory for FASTA file:", fout_fasta_dir, "\n")
        }, error = function(e) {
            stop(paste("Failed to create FASTA output directory:", fout_fasta_dir, "Error:", e$message))
        })
    }
    if (!dir.exists(fout_taxonomy_dir)) {
        tryCatch({
            dir.create(fout_taxonomy_dir, recursive = TRUE, showWarnings = FALSE)
            cat("Created output directory for taxonomy file:", fout_taxonomy_dir, "\n")
        }, error = function(e) {
            stop(paste("Failed to create taxonomy output directory:", fout_taxonomy_dir, "Error:", e$message))
        })
    }

    # Helper function for matching genera (adapted from dada2 R package)
    # Checks if gen.tax contains gen.binom, allowing for common delimiters or prefixes.
    matchGenera <- function(gen.tax, gen.binom, split.glyph="/") {
        if(is.na(gen.tax) || is.na(gen.binom)) { return(FALSE) }
        if((gen.tax == gen.binom) ||
           grepl(paste0("^", gen.binom, "[ _", split.glyph, "]"), gen.tax) ||
           grepl(paste0(split.glyph, gen.binom, "$"), gen.tax)) {
            return(TRUE)
        } else {
            return(FALSE)
        }
    }

    # Read and Parse Input FASTA ---
    xset <- Biostrings::DNAStringSet(Biostrings::readDNAStringSet(fin, format = "fasta"))

    if (length(xset) == 0) {
        warning("No sequences found in input. Empty output files will be created.")
        # Create empty but valid output files
        if (compress) {
            Biostrings::writeXStringSet(Biostrings::DNAStringSet(), fout_fasta, format="fasta", compress=TRUE)
        } else {
            Biostrings::writeXStringSet(Biostrings::DNAStringSet(), fout_fasta, format="fasta")
        }
        cat("ReferenceID\tTaxonomy\n", file = fout_taxonomy)
        return(invisible(TRUE))
    }

    descriptions <- names(xset)

    # Expected UNITE format: >Species_name|accession|reference_id|...|k__Fungi;p__...
    # Extracts the reference_id (typically the third field)
    ref_ids <- sapply(strsplit(descriptions, "\\|"), `[`, 3)

    if (any(is.na(ref_ids) | ref_ids == "")) {
        stop("FATAL: Failed to extract reference IDs from some sequence headers. Ensure they follow the UNITE format '...|reference_id|...'.")
    }

    # Handle duplicated sequence IDs by keeping only the first occurrence
    unique_ref_ids_logical <- !duplicated(ref_ids)
    if (any(!unique_ref_ids_logical)) {
        n_duplicates_removed <- sum(!unique_ref_ids_logical)
        warning(paste("WARNING:", n_duplicates_removed, "duplicated sequence IDs detected and removed (keeping first occurrence)."))
        xset <- xset[unique_ref_ids_logical]
        descriptions <- descriptions[unique_ref_ids_logical]
        ref_ids <- ref_ids[unique_ref_ids_logical]
    }
    names(xset) <- ref_ids

    # Extract and parse taxonomy strings (the last field in the description)
    taxl_parsed <- sapply(strsplit(descriptions, "\\|"), function(x) x[length(x)])
    taxa_list_from_fasta <- strsplit(taxl_parsed, ";")
    names(taxa_list_from_fasta) <- ref_ids

    # Construct and Clean Taxonomy Matrix ---
    num_ranks_unite <- 7 # Kingdom, Phylum, Class, Order, Family, Genus, Species
    unite_colnames <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    taxa.unite.mat <- matrix(NA_character_, nrow = length(ref_ids), ncol = num_ranks_unite)
    colnames(taxa.unite.mat) <- unite_colnames
    rownames(taxa.unite.mat) <- ref_ids

    # Populate matrix with taxonomic ranks, removing UNITE's 'r__' prefixes
    for (i in seq_along(taxa_list_from_fasta)) {
        len <- length(taxa_list_from_fasta[[i]])
        parsed_ranks <- sub("^[kpcofgs]__", "", taxa_list_from_fasta[[i]])
        taxa.unite.mat[i, 1:min(num_ranks_unite, len)] <- parsed_ranks[1:min(num_ranks_unite, len)]
    }

    # General cleaning of ambiguous terms across all ranks (excluding "Incertae_sedis")
    problematic_patterns <- "uncultured|unknown|unidentified|metagenome"
    taxa.unite.mat[grepl(problematic_patterns, taxa.unite.mat, ignore.case = TRUE)] <- NA_character_

    # Also clean empty strings or "s__" (empty species prefix) which may result from parsing
    taxa.unite.mat[taxa.unite.mat == "" | taxa.unite.mat == "s__"] <- NA_character_

    # Clean trailing "Incertae Sedis" terms specifically from Genus level upwards.
    # This loop iterates from right-to-left (Genus, Family, Order, etc.),
    # setting Incertae Sedis terms to NA until a valid non-Incertae Sedis rank is found.
    # The 'Species' column (column 7) is intentionally excluded from this specific cleanup
    # and will be handled by the 'Species Name Validation' block.
    genus_col_idx <- which(colnames(taxa.unite.mat) == "Genus") # Column index for Genus (expected to be 6)

    for (r_idx in 1:nrow(taxa.unite.mat)) {
        # Determine the actual rightmost non-NA column in the current row
        actual_rightmost_non_na <- max(which(!is.na(taxa.unite.mat[r_idx, ])), 0)

        # Set the starting column index for backward iteration for this cleanup.
        # It's the minimum of the actual rightmost non-NA rank's index and the Genus column index (6).
        # This ensures we always start at or above Genus, excluding Species.
        start_c_idx <- min(actual_rightmost_non_na, genus_col_idx)

        # Skip row if there are no non-NA ranks up to Genus
        if (start_c_idx == 0) next

        # Iterate backwards from the determined 'start_c_idx' down to Kingdom (column 1)
        for (c_idx in start_c_idx:1) {
            # If the current rank is not NA and contains "Incertae_sedis" (or variations)
            if (!is.na(taxa.unite.mat[r_idx, c_idx]) &&
                grepl("Incertae_sedis.*", taxa.unite.mat[r_idx, c_idx], ignore.case = TRUE)) {
                taxa.unite.mat[r_idx, c_idx] <- NA_character_ # Set to NA
            } else {
                # Stop cleaning for this row when a valid, non-"Incertae Sedis" rank is found
                break
            }
        }
    }

    # Propagate NAs to maintain strict taxonomic hierarchy.
    # If a higher rank (e.g., Phylum) is NA, all subsequent lower ranks (Class, Order, etc.) are also set to NA.
    for (c_idx in 1:(ncol(taxa.unite.mat) - 1)) {
        na_rows <- is.na(taxa.unite.mat[, c_idx])
        taxa.unite.mat[na_rows, (c_idx + 1):ncol(taxa.unite.mat)] <- NA_character_
    }

    # Species Name Validation ---
    # This block specifically validates and processes the Species rank.
    # It checks for genus-species consistency, generic species terms ("sp"),
    # and other ambiguous terms within the unite file species string.

    if (include.species) {
        # Extract raw species field, removing 's__' prefix
        raw_species_field <- sapply(taxa_list_from_fasta, function(x) if(length(x) >= 7) sub("s__", "", x[7]) else NA_character_)
        # Get the cleaned genus name from the matrix (reflecting prior cleanups)
        genus_from_col6 <- taxa.unite.mat[, "Genus"]

        # Clean genus and raw species fields for matching and term identification
        genus_clean_for_match <- gsub("Candidatus |\\[|\\]", "", genus_from_col6)
        binom_raw_clean <- gsub("Candidatus |\\[|\\]", "", raw_species_field)

        # Split the raw species string into potential genus and epithet parts
        # Tries to split by underscore first, then by space.
        binom_parsed_parts_list <- lapply(binom_raw_clean, function(s_name) {
            if (is.na(s_name)) return(c(NA_character_, NA_character_))
            parts <- if (grepl("_", s_name)) unlist(strsplit(s_name, "_")) else unlist(strsplit(s_name, " "))
            c(parts[1], if(length(parts) >= 2) parts[2] else NA_character_)
        })
        binom_parsed_parts <- do.call(rbind, binom_parsed_parts_list)
        colnames(binom_parsed_parts) <- c("parsed_genus", "parsed_epithet")

        # Check if the parsed genus from the species name matches the genus in the taxonomy matrix
        gen.match <- rep(FALSE, nrow(taxa.unite.mat))
        can_match_idx <- !is.na(genus_clean_for_match) & !is.na(binom_parsed_parts[, "parsed_genus"])

        if(any(can_match_idx)){
            gen.match[can_match_idx] <- mapply(matchGenera,
                                               genus_clean_for_match[can_match_idx],
                                               binom_parsed_parts[can_match_idx, "parsed_genus"],
                                               MoreArgs = list(split.glyph = "-"))
        }

        # Check for NA epithet
        is.NA_epithet <- is.na(binom_parsed_parts[, "parsed_epithet"])

        # Identify generic species terms (e.g., "sp", "spp") - uses exact match for "sp"
        is.sp <- grepl("^sp$", binom_parsed_parts[, "parsed_epithet"], ignore.case = TRUE)
        is.sp[is.na(is.sp)] <- FALSE # Ensure no NA values in logical vector

        # Identify other problematic terms in the raw species string (e.g., uncultured, Incertae_sedis)
        is.problem_term <- grepl("uncultured|unidentified|endosymbiont|Incertae_sedis.*", binom_raw_clean, ignore.case = TRUE)
        is.problem_term[is.na(is.problem_term)] <- FALSE # Ensure no NA values in logical vector

        # Determine overall validity of the species name:
        # Requires genus match, non-NA epithet, not a generic 'sp' term, and not other problematic terms.
        valid.spec <- gen.match & !is.NA_epithet & !is.sp & !is.problem_term
        valid.spec[is.na(valid.spec)] <- FALSE # Crucial for cases where components of valid.spec might be NA

        # Construct the full species name (Genus_epithet) for valid entries, otherwise set to NA
        full_species_to_store <- rep(NA_character_, length(valid.spec))
        can_form <- valid.spec & !is.na(genus_from_col6) & !is.na(binom_parsed_parts[, "parsed_epithet"])
        if(any(can_form)) {
            full_species_to_store[can_form] <- paste0(genus_from_col6[can_form], "_", binom_parsed_parts[can_form, "parsed_epithet"])
        }
        taxa.unite.mat[, "Species"] <- full_species_to_store
    } else {
        # If include.species is FALSE, set the Species column to NA for all entries
        taxa.unite.mat[, "Species"] <- NA_character_
    }

    # Write Output Files ---

    # Filter sequences to keep only those with processed taxonomy
    seqs.keep.ids <- rownames(taxa.unite.mat)
    seqs.keep <- xset[seqs.keep.ids]

    if (length(seqs.keep) != length(seqs.keep.ids)) {
        stop("FATAL: Mismatch in the number of sequences and IDs after processing. Internal error.")
    }

    # Write fasta file (compressed if requested)
    if (length(seqs.keep) > 0) {
        if (compress) {
            Biostrings::writeXStringSet(seqs.keep, fout_fasta, format="fasta", width=20000L, compress=TRUE)
        } else {
            Biostrings::writeXStringSet(seqs.keep, fout_fasta, format="fasta", width=20000L)
        }
    } else {
        warning("No valid sequences remained after filtering. Empty FASTA file created.")
        cat(">", file = fout_fasta) # Ensure a valid empty FASTA header is written for empty file
    }

    # Prepare and write taxonomy TSV
    # Select columns based on 'include.species' flag
    output_tax_mat <- if (!include.species) taxa.unite.mat[, 1:6, drop = FALSE] else taxa.unite.mat
    prefixes <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")

    # Helper function to format a single row of taxonomy matrix into the desired string format
    format_taxonomy_row <- function(tax_row) {
        num_levels <- length(tax_row)
        current_prefixes <- prefixes[1:num_levels]
        parts <- character(num_levels)
        first_na_hit <- FALSE # Flag to ensure all subsequent ranks are also NA if one is NA
        for (i in 1:num_levels) {
            if (first_na_hit || is.na(tax_row[i])) {
                parts[i] <- current_prefixes[i] # If NA or a prior rank was NA, just output prefix
                first_na_hit <- TRUE
            } else {
                parts[i] <- paste0(current_prefixes[i], tax_row[i]) # Output prefix + taxonomy name
            }
        }
        paste(parts, collapse = ";")
    }

    if (nrow(output_tax_mat) > 0) {
        formatted_taxonomy <- apply(output_tax_mat, 1, format_taxonomy_row)
        taxonomy_df <- data.frame(ReferenceID = rownames(output_tax_mat), Taxonomy = formatted_taxonomy)
        write.table(taxonomy_df, fout_taxonomy, sep = "\t", quote = FALSE, row.names = FALSE)
    } else {
        warning("No valid taxonomy remained after filtering. Taxonomy file with only headers created.")
        cat("ReferenceID\tTaxonomy\n", file = fout_taxonomy) # Write header for empty file
    }

    # Final Summary ---

    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    cat("\nProcessing Summary (UNITE):\n")
    cat("----------------------------\n")
    cat(nrow(taxa.unite.mat), "reference sequences were selected for output.\n")
    cat("Breakdown by Kingdom:\n")
    print(table(taxa.unite.mat[, "Kingdom"], useNA = "ifany"))
    if (include.species) {
        species_count <- sum(!is.na(taxa.unite.mat[, "Species"]))
        cat(species_count, "entries include validated species identity.\n")
    }
    cat("Execution time:", round(execution_time, 1), "seconds.\n")
    cat("----------------------------\n\n")

    return(invisible(TRUE))
}
