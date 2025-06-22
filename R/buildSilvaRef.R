# Biostrings required for fasta file handling
library(Biostrings)

#' Create sequence and taxonomy file for microbiome analysis
#'
#' Reads a SILVA fasta (tested on NR SSU set) and corresponding taxonomic map file, extracts
#' annotations, filters them, and outputs a refined fasta and taxonomy file. Optionally
#' includes a set number of eukaryotic sequences. The matchGenera function is adapted
#' from the dada2 R package.
#'
#' @importFrom Biostrings DNAStringSet readDNAStringSet readRNAStringSet writeXStringSet
#' @importFrom utils read.table write.table
#'
#' @param fin Path to the input SILVA fasta file (RNA sequences).
#' @param ftax Path to the SILVA taxonomy map file with valid taxon strings and levels.
#' @param fout_fasta Path to the output fasta file.
#' @param fout_taxonomy Path to the output taxonomy TSV file.
#' @param include.species Logical; whether to include species names in taxonomy.
#' @param compress Logical; whether to compress output sequences using gzip.
#' @param n_euk Integer; number of random eukaryotic sequences to include in the final output.
#'
#' @return Writes a fasta and a taxonomy file to disk. Returns nothing.
#' @export
#'
#' @examples
#' \dontrun{
#' buildSilvaRef("SILVA.fasta", "SILVA_taxonomy.tsv",
#'                      "output.fasta", "output_taxonomy.tsv",
#'                      include.species = TRUE, compress = FALSE,
#'                      n_euk = 100)
#' }

buildSilvaRef <- function(fin, ftax, fout_fasta, fout_taxonomy,
                          include.species = TRUE,
                          compress = FALSE,
                          n_euk = 100) {

    # Start timing
    start_time <- Sys.time()

    # Parameter validation
    if (missing(fin) || !is.character(fin) || length(fin) != 1) {
        stop("'fin' must be a single character string specifying the input fasta file path")
    }

    if (missing(ftax) || !is.character(ftax) || length(ftax) != 1) {
        stop("'ftax' must be a single character string specifying the taxonomy map file path")
    }

    if (missing(fout_fasta) || !is.character(fout_fasta) || length(fout_fasta) != 1) {
        stop("'fout_fasta' must be a single character string specifying the output fasta file path")
    }

    if (missing(fout_taxonomy) || !is.character(fout_taxonomy) || length(fout_taxonomy) != 1) {
        stop("'fout_taxonomy' must be a single character string specifying the output taxonomy file path")
    }

    # Check input file existence
    if (!file.exists(fin)) {
        stop(paste("Input fasta file does not exist:", fin))
    }
    if (!file.exists(ftax)) {
        stop(paste("Taxonomy map file does not exist:", ftax))
    }

    # Create output directories if needed
    fout_fasta_dir <- dirname(fout_fasta)
    fout_taxonomy_dir <- dirname(fout_taxonomy)

    # Create fasta output directory
    if (!dir.exists(fout_fasta_dir)) {
        cat("Creating output directory for fasta file:", fout_fasta_dir, "\n")
        tryCatch({
            dir.create(fout_fasta_dir, recursive = TRUE, showWarnings = FALSE)
        }, error = function(e) {
            stop(paste("Failed to create output directory for fasta file:", fout_fasta_dir,
                       "Error:", e$message))
        })
    }

    # Create taxonomy output directory
    if (!dir.exists(fout_taxonomy_dir)) {
        cat("Creating output directory for taxonomy file:", fout_taxonomy_dir, "\n")
        tryCatch({
            dir.create(fout_taxonomy_dir, recursive = TRUE, showWarnings = FALSE)
        }, error = function(e) {
            stop(paste("Failed to create output directory for taxonomy file:", fout_taxonomy_dir,
                       "Error:", e$message))
        })
    }

    # Validate logical parameters
    if (!is.logical(include.species) || length(include.species) != 1) {
        stop("'include.species' must be a single logical value (TRUE or FALSE)")
    }
    if (!is.logical(compress) || length(compress) != 1) {
        stop("'compress' must be a single logical value (TRUE or FALSE)")
    }

    # Validate n_euk parameter
    if (!is.numeric(n_euk) || length(n_euk) != 1 || n_euk < 0 || n_euk != as.integer(n_euk)) {
        stop("'n_euk' must be a single non-negative integer")
    }

    # Helper function for matching genera. The code adapted from dada2 R package
    matchGenera <- function(gen.tax, gen.binom, split.glyph="/") {
        if(is.na(gen.tax) || is.na(gen.binom)) {
            return(FALSE)
        }
        if((gen.tax==gen.binom) ||
           grepl(paste0("^", gen.binom, "[ _", split.glyph, "]"), gen.tax) ||
           grepl(paste0(split.glyph, gen.binom, "$"), gen.tax)) {
            return(TRUE)
        } else {
            return(FALSE)
        }
    }

    # Read and convert sequences
    xset <- DNAStringSet(readRNAStringSet(fin, format = "fasta"))

    # Extract sequence metadata
    descriptions <- names(xset)
    ref_ids <- sapply(strsplit(descriptions, "\\s"), `[` , 1)

    # Handle duplicate sequence IDs
    unique_ref_ids_logical <- !duplicated(ref_ids)
    if (any(!unique_ref_ids_logical)) {
        n_duplicates_removed <- sum(!unique_ref_ids_logical)
        warning(paste("WARNING:", n_duplicates_removed, "duplicated sequence IDs detected and removed (keeping first occurrence)."))

        xset <- xset[unique_ref_ids_logical]
        descriptions <- descriptions[unique_ref_ids_logical]
        ref_ids <- sapply(strsplit(descriptions, "\\s"), `[` , 1)

        if (any(duplicated(ref_ids))) {
            stop("FATAL: Duplicated sequence IDs still detected after filtering. Please check input fasta")
        }
    }

    names(xset) <- ref_ids

    # Parse taxonomic annotations from fasta
    taxl_parsed <- gsub("^[A-Za-z0-9.]+\\s", "", descriptions)
    taxl_parsed <- gsub(";YM;", ";", taxl_parsed)
    taxa_list_from_fasta <- strsplit(taxl_parsed, ";")
    names(taxa_list_from_fasta) <- ref_ids

    # Read SILVA taxonomy reference
    silva.taxa <- read.table(ftax, sep = "\t",
                             col.names = c("Taxon", "V2", "Level", "V4", "V5"),
                             stringsAsFactors = FALSE)[, c("Taxon", "Level")]

    # Filter by kingdom
    kingdom_from_fasta <- sapply(taxa_list_from_fasta, function(t) if(length(t) > 0) t[1] else NA_character_)

    is_ba_from_fasta <- kingdom_from_fasta %in% c("Bacteria", "Archaea")
    is_ba_from_fasta[is.na(is_ba_from_fasta)] <- FALSE

    ref_ids.ba <- names(taxa_list_from_fasta[is_ba_from_fasta])

    # Construct Bacteria/Archaea taxonomy matrix
    num_ranks_ba <- 6
    ba_colnames <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
    taxa.ba.mat <- matrix(NA_character_, nrow = length(ref_ids.ba), ncol = num_ranks_ba)
    colnames(taxa.ba.mat) <- ba_colnames

    if (length(ref_ids.ba) > 0) {
        rownames(taxa.ba.mat) <- ref_ids.ba
        current_taxa_ba_list <- taxa_list_from_fasta[ref_ids.ba]

        for (i in seq_along(current_taxa_ba_list)) {
            len <- length(current_taxa_ba_list[[i]])
            taxa.ba.mat[i, 1:min(num_ranks_ba, len)] <- current_taxa_ba_list[[i]][1:min(num_ranks_ba, len)]
        }
    }

    # Validate taxonomy against SILVA reference
    if (nrow(taxa.ba.mat) > 0) {
        taxa.ba.mat.string <- matrix(NA_character_, nrow = nrow(taxa.ba.mat), ncol = ncol(taxa.ba.mat))
        for (r in 1:nrow(taxa.ba.mat)) {
            current_path <- ""
            for (c in 1:ncol(taxa.ba.mat)) {
                if (!is.na(taxa.ba.mat[r,c]) && taxa.ba.mat[r,c] != "") {
                    current_path <- paste0(current_path, taxa.ba.mat[r,c], ";")
                    taxa.ba.mat.string[r,c] <- current_path
                } else {
                    taxa.ba.mat[r,c] <- NA_character_
                    break
                }
            }
        }

        taxa.ba.mat.is_valid <- matrix(TRUE, nrow = nrow(taxa.ba.mat.string), ncol = ncol(taxa.ba.mat.string))
        for (r_idx in 1:nrow(taxa.ba.mat.string)) {
            for (c_idx in 1:ncol(taxa.ba.mat.string)) {
                current_path_val <- taxa.ba.mat.string[r_idx, c_idx]
                if (!is.na(current_path_val)) {
                    if (!(current_path_val %in% silva.taxa$Taxon)) {
                        taxa.ba.mat.is_valid[r_idx, c_idx] <- FALSE
                    }
                }
            }
        }

        for (r_idx in 1:nrow(taxa.ba.mat)) {
            for (c_idx in 1:ncol(taxa.ba.mat)) {
                if (!taxa.ba.mat.is_valid[r_idx, c_idx]) {
                    taxa.ba.mat[r_idx, c_idx:ncol(taxa.ba.mat)] <- NA_character_
                    break
                }
            }
        }

        # Clean taxonomy terms
        taxa.ba.mat[grepl("uncultured|unknown|unidentified", taxa.ba.mat, ignore.case = TRUE)] <- NA_character_

        # Clean trailing Incertae Sedis terms
        for (r_idx in 1:nrow(taxa.ba.mat)) {
            rightmost_non_na <- max(which(!is.na(taxa.ba.mat[r_idx, ])), 0)
            if (rightmost_non_na > 0) {
                for (c_idx in rightmost_non_na:1) {
                    if (!is.na(taxa.ba.mat[r_idx, c_idx]) &&
                        grepl("Incertae Sedis", taxa.ba.mat[r_idx, c_idx], ignore.case = FALSE)) {
                        taxa.ba.mat[r_idx, c_idx] <- NA_character_
                    } else {
                        break
                    }
                }
            }
        }

        for (r_idx in 1:nrow(taxa.ba.mat)) {
            for (c_idx in 1:(ncol(taxa.ba.mat)-1)) {
                if(is.na(taxa.ba.mat[r_idx, c_idx])) {
                    taxa.ba.mat[r_idx, (c_idx+1):ncol(taxa.ba.mat)] <- NA_character_
                }
            }
        }
    }

    # Handle species information
    species_col_name <- "Species"
    final_colnames <- ba_colnames

    if (include.species) {
        final_colnames <- c(ba_colnames, species_col_name)

        species_values_col <- matrix(NA_character_, nrow = nrow(taxa.ba.mat), ncol = 1)
        colnames(species_values_col) <- species_col_name
        if(nrow(taxa.ba.mat) > 0) rownames(species_values_col) <- rownames(taxa.ba.mat)
        taxa.ba.mat <- cbind(taxa.ba.mat, species_values_col)

        if (nrow(taxa.ba.mat) > 0) {
            raw_species_field <- sapply(taxa_list_from_fasta[rownames(taxa.ba.mat)],
                                        function(x) if(length(x) >= 7) x[7] else NA_character_)

            genus_from_col6 <- taxa.ba.mat[, "Genus"]
            genus_clean_for_match <- gsub("Candidatus ", "", genus_from_col6)
            genus_clean_for_match <- gsub("\\[|\\]", "", genus_clean_for_match)

            binom_raw_clean <- gsub("Candidatus ", "", raw_species_field)
            binom_raw_clean <- gsub("\\[|\\]", "", binom_raw_clean)

            binom_parsed_parts_list <- lapply(strsplit(binom_raw_clean, "\\s"), function(parts) {
                c(if(length(parts) >= 1) parts[1] else NA_character_,
                  if(length(parts) >= 2) parts[2] else NA_character_)
            })
            binom_parsed_parts <- do.call(rbind, binom_parsed_parts_list)
            colnames(binom_parsed_parts) <- c("parsed_genus", "parsed_epithet")

            gen.match <- rep(FALSE, nrow(taxa.ba.mat))
            can_match_idx <- !is.na(genus_clean_for_match) & genus_clean_for_match != "" &
                !is.na(binom_parsed_parts[, "parsed_genus"]) & binom_parsed_parts[, "parsed_genus"] != ""

            if(any(can_match_idx)){
                gen.match[can_match_idx] <- mapply(matchGenera,
                                                   genus_clean_for_match[can_match_idx],
                                                   binom_parsed_parts[can_match_idx, "parsed_genus"],
                                                   MoreArgs = list(split.glyph = "-"))
            }

            is.NA_epithet <- is.na(binom_parsed_parts[, "parsed_epithet"]) | binom_parsed_parts[, "parsed_epithet"] == ""
            is.sp <- grepl("sp\\.", binom_parsed_parts[, "parsed_epithet"])
            is.sp[is.na(is.sp)] <- FALSE

            is.endo <- (grepl("endosymbiont", binom_parsed_parts[, "parsed_genus"], ignore.case = TRUE) |
                            grepl("endosymbiont", binom_parsed_parts[, "parsed_epithet"], ignore.case = TRUE))
            is.endo[is.na(is.endo)] <- FALSE

            is.uncult <- (grepl("[Uu]ncultured", binom_parsed_parts[, "parsed_genus"]) |
                              grepl("[Uu]ncultured", binom_parsed_parts[, "parsed_epithet"]))
            is.uncult[is.na(is.uncult)] <- FALSE

            is.unident <- (grepl("[Uu]nidentified", binom_parsed_parts[, "parsed_genus"]) |
                               grepl("[Uu]nidentified", binom_parsed_parts[, "parsed_epithet"]))
            is.unident[is.na(is.unident)] <- FALSE

            valid.spec <- gen.match & !is.NA_epithet & !is.sp & !is.endo & !is.uncult & !is.unident
            valid.spec[is.na(valid.spec)] <- FALSE

            full_species_to_store <- character(length(valid.spec))
            full_species_to_store[] <- NA_character_

            can_form_full_species <- valid.spec &
                !is.na(genus_from_col6) & genus_from_col6 != "" &
                !is.na(binom_parsed_parts[, "parsed_epithet"]) & binom_parsed_parts[, "parsed_epithet"] != ""

            if(any(can_form_full_species)) {
                full_species_to_store[can_form_full_species] <- paste0(
                    genus_from_col6[can_form_full_species], "_",
                    binom_parsed_parts[can_form_full_species, "parsed_epithet"]
                )
            }

            taxa.ba.mat[, species_col_name] <- full_species_to_store
        }
    }

    # Sample eukaryotic sequences
    euk_indices <- which(kingdom_from_fasta == "Eukaryota" & !is.na(kingdom_from_fasta))
    euk_ids_available <- ref_ids[euk_indices]
    n_euk_actual <- 0
    euk.keep_ids <- character(0)

    if (n_euk > 0 && length(euk_ids_available) > 0) {
        if (length(euk_ids_available) < n_euk) {
            warning(paste("Requested", n_euk, "eukaryotic sequences, but only", length(euk_ids_available), "are available. Using all available."))
            n_euk_actual <- length(euk_ids_available)
        } else {
            n_euk_actual <- n_euk
        }
        set.seed(100)
        euk.keep_ids <- sample(euk_ids_available, n_euk_actual)
    }

    taxa.euk.mat <- matrix(NA_character_, nrow = n_euk_actual, ncol = length(final_colnames))
    colnames(taxa.euk.mat) <- final_colnames
    if (n_euk_actual > 0) {
        rownames(taxa.euk.mat) <- euk.keep_ids
        taxa.euk.mat[, "Kingdom"] <- "Eukaryota"
    }

    # Combine taxonomy matrices
    taxa.mat.final <- rbind(taxa.ba.mat, taxa.euk.mat)

    # Prepare sequences for output
    seqs.keep.ids <- rownames(taxa.mat.final)

    if(length(seqs.keep.ids) > 0) {
        if (any(is.na(seqs.keep.ids))) stop("FATAL: NA found in final sequence IDs before output.")
        if (!all(seqs.keep.ids %in% ref_ids)) stop("FATAL: Some final sequence IDs do not originate from input ref_ids.")
    }

    seqs.keep <- xset[seqs.keep.ids]
    if(length(seqs.keep) != length(seqs.keep.ids)) stop("FATAL: Mismatch in number of sequences to keep and their IDs.")

    # Write FASTA output
    if (length(seqs.keep) > 0) {
        if (compress) {
            writeXStringSet(seqs.keep, fout_fasta, format = "fasta", width = 20000L, compress = TRUE)
        } else {
            writeXStringSet(seqs.keep, fout_fasta, format = "fasta", width = 20000L)
        }
    } else {
        cat(">", file=fout_fasta)
        warning("No sequences selected for FASTA output. Empty FASTA file created.")
    }

    # Prepare and write taxonomy TSV
    if (nrow(taxa.mat.final) > 0) {
        prefixes <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")

        # Define the taxonomy formatting function separately
        format_taxonomy_row <- function(tax_row) {
            num_levels <- length(tax_row)
            current_prefixes <- prefixes[1:num_levels]
            parts <- character(num_levels)

            first_na_or_empty_hit <- FALSE
            for (i in 1:num_levels) {
                if (first_na_or_empty_hit) {
                    parts[i] <- current_prefixes[i]
                } else if (!is.na(tax_row[i]) && tax_row[i] != "") {
                    parts[i] <- paste0(current_prefixes[i], tax_row[i])
                } else {
                    parts[i] <- current_prefixes[i]
                    first_na_or_empty_hit <- TRUE
                }
            }
            paste(parts, collapse = ";")
        }

        formatted_taxonomy_strings <- apply(taxa.mat.final, 1, format_taxonomy_row)

        taxonomy_output_df <- data.frame(
            ReferenceID = rownames(taxa.mat.final),
            Taxonomy = formatted_taxonomy_strings,
            stringsAsFactors = FALSE
        )

        write.table(taxonomy_output_df, file = fout_taxonomy,
                    sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    } else {
        cat("ReferenceID\tTaxonomy\n", file = fout_taxonomy)
        warning("No sequences selected for taxonomy output. Taxonomy file with only headers created.")
    }

    # Calculate execution time and report summary
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    cat("\nProcessing Summary:\n")
    cat("--------------------\n")
    if(nrow(taxa.mat.final) > 0){
        cat(nrow(taxa.mat.final), "reference sequences were selected for output.\n")
        cat("Breakdown by Kingdom (from final processed taxonomy):\n")
        if("Kingdom" %in% colnames(taxa.mat.final)){
            print(table(taxa.mat.final[, "Kingdom"], useNA = "ifany"))
        } else {
            cat("Kingdom column not found in final taxonomy matrix for breakdown.\n")
        }
        if (include.species && "Species" %in% colnames(taxa.mat.final)) {
            species_count <- sum(!is.na(taxa.mat.final[, "Species"]) & taxa.mat.final[, "Species"] != "")
            cat(species_count, "entries include species identity\n")
        }
    } else {
        cat("No sequences were selected for output.\n")
    }
    cat("Execution time:", round(execution_time, 1), "seconds\n")
    cat("--------------------\n\n")

    return(invisible(TRUE))
}
