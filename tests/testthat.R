# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(taxreftool)
test_check("taxreftool")


test_that("buildUniteRef creates expected files", {
    # Create a dummy fasta file for testing
    temp_fasta <- tempfile(fileext = ".fasta")
    temp_out_fasta <- tempfile(fileext = ".fasta")
    temp_out_taxonomy <- tempfile(fileext = ".tsv")

    # Write a minimal dummy FASTA content
    # This needs to be a valid UNITE-like header for your function to parse it
    dummy_content <- paste0(
        ">SH100.00FU|unite_1|SH100.00FU|k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__;g__;s__\n",
        "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT\n",
        ">SH101.00FU|unite_2|SH101.00FU|k__Fungi;p__Basidiomycota;c__;o__;f__;g__;s__\n",
        "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT\n"
    )
    writeLines(dummy_content, temp_fasta)

    # Run your function
    buildUniteRef(fin = temp_fasta,
                  fout_fasta = temp_out_fasta,
                  fout_taxonomy = temp_out_taxonomy,
                  include.species = TRUE,
                  compress = FALSE)

    # Assertions
    expect_true(file.exists(temp_out_fasta))
    expect_true(file.exists(temp_out_taxonomy))

    # Check if output FASTA is not empty (and perhaps check contents)
    output_fasta_seqs <- Biostrings::readDNAStringSet(temp_out_fasta)
    expect_gt(length(output_fasta_seqs), 0)

    # Check taxonomy file contents (e.g., number of rows)
    output_tax <- read.delim(temp_out_taxonomy, stringsAsFactors = FALSE)
    expect_equal(nrow(output_tax), 2)
    expect_true(all(c("ReferenceID", "Taxonomy") %in% colnames(output_tax)))
    expect_equal(output_tax$ReferenceID[1], "SH100.00FU")
    # Add more specific checks for parsed taxonomy as needed

    # Clean up
    unlink(c(temp_fasta, temp_out_fasta, temp_out_taxonomy))
})

test_that("buildUniteRef handles empty input", {
    temp_fasta_empty <- tempfile(fileext = ".fasta")
    file.create(temp_fasta_empty) # Create an empty file
    temp_out_fasta <- tempfile(fileext = ".fasta")
    temp_out_taxonomy <- tempfile(fileext = ".tsv")

    # Expect a warning for empty input
    expect_warning(buildUniteRef(fin = temp_fasta_empty,
                                 fout_fasta = temp_out_fasta,
                                 fout_taxonomy = temp_out_taxonomy))

    expect_true(file.exists(temp_out_fasta))
    expect_true(file.exists(temp_out_taxonomy))
    expect_equal(length(Biostrings::readDNAStringSet(temp_out_fasta)), 0)

    unlink(c(temp_fasta_empty, temp_out_fasta, temp_out_taxonomy))
})
