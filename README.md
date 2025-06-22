# taxreftool: Tools for Building Taxonomic Reference Databases

[![R package](https://img.shields.io/badge/R-package-blue.svg)](https://www.r-project.org/) [![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

## Overview

The `taxreftool` package provides utilities for processing and formatting taxonomic reference databases from SILVA and UNITE releases. Tools are aimed for microbial community analysis workflows, particularly for amplicon sequencing studies targeting bacterial, archaeal, and fungal communities.

The package implements taxonomic validation procedures, including:

-   Hierarchical taxonomy validation against reference standards

-   Species name validation using genus-epithet matching algorithms

-   Systematic removal of ambiguous taxonomic annotations

-   Quality control measures for sequence and metadata consistency

## Key Features

-   **SILVA Processing**: Comprehensive processing of SILVA SSU rRNA reference sequences with taxonomic validation

-   **UNITE Processing**: Specialized handling of UNITE fungal ITS reference sequences

-   **Taxonomic Validation**: Implementation of genus-species matching algorithms adapted from the dada2 package

-   **Quality Control**: Systematic removal of problematic annotations (uncultured, unidentified, Incertae sedis)

-   **Flexible Output**: Support for compressed outputs and configurable species inclusion

## Data Sources

### SILVA Database Files

The SILVA SSU reference sequences and taxonomy files can be obtained from the SILVA database project website ([https://www.arb-silva.de/archive/).](https://www.arb-silva.de/download/).) Download the "SILVA SSU Ref NR" FASTA file (e.g., `SILVA_138.2_SSURef_NR99_tax_silva.fasta`) and the corresponding taxonomy map file (e.g., `tax_slv_ssu_138.2.txt`) from the same release version. Both files are required for processing with `buildSilvaRef()`.

### UNITE Database Files

The UNITE fungal ITS reference sequences are available from the UNITE database (<https://unite.ut.ee/repository.php>) Download the "General FASTA release" file (e.g., `sh_general_release_dynamic_25.07.2023.fasta`) which contains both sequence data and embedded taxonomic information. The UNITE format includes taxonomy strings within the FASTA headers, eliminating the need for separate taxonomy files.

## Installation

### Prerequisites

The package requires R (≥ 4.0.0) and the following dependencies:

```         
# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")
```

### Package Installation

#### From Source (Recommended)

```         
# Install devtools if not already installed
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

# Install from local source
devtools::install_local("path/to/taxreftool", dependencies = TRUE)
```

#### Manual Installation

```         
# From command line
R CMD INSTALL taxreftool_1.0.0.tar.gz
```

## Functions

### `buildSilvaRef()`

Processes SILVA reference sequences for bacterial and archaeal taxonomy assignment.

#### Methodology

The function implements a multi-step validation process:

1.  **Sequence Processing**: Converts RNA sequences to DNA format and handles duplicate sequence identifiers

2.  **Taxonomic Parsing**: Extracts taxonomic information from FASTA headers and validates against SILVA taxonomy maps

3.  **Hierarchy Validation**: Ensures taxonomic consistency by validating each rank against the SILVA reference taxonomy

4.  **Quality Filtering**: Removes sequences with problematic annotations (uncultured, unidentified, Incertae sedis)

5.  **Species Validation**: Implements genus-epithet matching algorithms for species-level assignments

6.  **Eukaryotic Sampling**: Optional inclusion of random eukaryotic sequences for comprehensive databases

#### Usage

```         
library(taxreftool)

# Basic usage
buildSilvaRef(
    fin = "SILVA_138_SSURef_NR99_tax_silva.fasta",
    ftax = "tax_slv_ssu_138.txt",
    fout_fasta = "silva_formatted.fasta",
    fout_taxonomy = "silva_taxonomy.tsv"
)

# Advanced usage with options
buildSilvaRef(
    fin = "SILVA_138_SSURef_NR99_tax_silva.fasta",
    ftax = "tax_slv_ssu_138.txt", 
    fout_fasta = "silva_formatted.fasta.gz",
    fout_taxonomy = "silva_taxonomy.tsv",
    include.species = TRUE,
    compress = TRUE,
    n_euk = 500
)
```

#### Parameters

-   `fin`: Path to input SILVA FASTA file (RNA sequences)

-   `ftax`: Path to SILVA taxonomy map file

-   `fout_fasta`: Output FASTA file path

-   `fout_taxonomy`: Output taxonomy TSV file path

-   `include.species`: Include species-level assignments (default: TRUE)

-   `compress`: Compress output FASTA with gzip (default: FALSE)

-   `n_euk`: Number of eukaryotic sequences to include (default: 100)

### `buildUniteRef()`

Processes UNITE reference sequences for fungal taxonomy assignment.

#### Methodology

The function processes UNITE releases through:

1.  **Header Parsing**: Extracts reference IDs and taxonomic strings from UNITE-formatted headers

2.  **Duplicate Handling**: Identifies and removes duplicate sequence entries

3.  **Taxonomic Matrix Construction**: Builds hierarchical taxonomy matrices with rank-specific validation

4.  **Incertae Sedis Cleaning**: Systematic removal of trailing "Incertae sedis" annotations while preserving valid taxonomy

5.  **Species Validation**: Validates species names through genus-epithet consistency checking

6.  **Output Formatting**: Generates standardized taxonomy strings with QIIME-compatible prefixes

#### Usage

```         
library(taxreftool)

# Basic usage
buildUniteRef(
    fin = "sh_general_release_dynamic_25.07.2023.fasta",
    fout_fasta = "unite_formatted.fasta",
    fout_taxonomy = "unite_taxonomy.tsv"
)

# With compression and species validation
buildUniteRef(
    fin = "sh_general_release_dynamic_25.07.2023.fasta",
    fout_fasta = "unite_formatted.fasta.gz",
    fout_taxonomy = "unite_taxonomy.tsv",
    include.species = TRUE,
    compress = TRUE
)
```

#### Parameters

-   `fin`: Path to input UNITE FASTA file

-   `fout_fasta`: Output FASTA file path

-   `fout_taxonomy`: Output taxonomy TSV file path

-   `include.species`: Include validated species names (default: TRUE)

-   `compress`: Compress output FASTA with gzip (default: FALSE)

## Output Format

Both functions generate two output files:

### FASTA File

-   Standard FASTA format with reference IDs as headers

-   Optional gzip compression

-   Sequences ready for taxonomic assignment workflows

### Taxonomy File

-   Tab-separated format with columns: `ReferenceID` and `Taxonomy`

-   QIIME-compatible taxonomy strings with rank prefixes:

    -   `k__` (Kingdom), `p__` (Phylum), `c__` (Class), `o__` (Order), `f__` (Family), `g__` (Genus), `s__` (Species)

-   Hierarchical consistency maintained throughout

Example taxonomy string:

```         
k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Hypocreales;f__Nectriaceae;g__Fusarium;s__Fusarium_oxysporum
```

## Quality Control Features

### Taxonomic Validation

-   Cross-reference validation against authoritative taxonomy databases

-   Hierarchical consistency checking

-   Systematic removal of ambiguous terms

### Species Name Validation

-   Genus-epithet matching using algorithms adapted from dada2

-   Identification and removal of generic species terms (e.g., "sp.")

-   Validation of binomial nomenclature consistency

### Problematic Annotation Handling

-   Systematic removal of terms: "uncultured", "unidentified", "unknown", "metagenome"

-   Intelligent handling of "Incertae sedis" annotations

-   Preservation of valid taxonomy while removing ambiguous classifications

## Performance Considerations

-   **Memory Usage**: Functions are optimized for large reference databases (\>1M sequences)

-   **Processing Time**: Typical processing times range from 30 seconds to several minutes depending on database size

-   **Disk Space**: Consider output compression for large databases

## Integration with Analysis Pipelines

The formatted outputs are compatible with major microbial community analysis tools:

-   **QIIME 2**: Direct import of taxonomy files for feature classification

-   **dada2**: Compatible reference databases for `assignTaxonomy()` and `assignSpecies()`

-   **mothur**: Formatted databases for taxonomic classification workflows

-   **Minimap2**: Reference databases for sequence alignment based classification

## Example Workflow

```         
library(taxreftool)

# Process SILVA database for bacteria/archaea
buildSilvaRef(
    fin = "SILVA_138_SSURef_NR99_tax_silva.fasta",
    ftax = "tax_slv_ssu_138.txt",
    fout_fasta = "references/silva_138_formatted.fasta.gz",
    fout_taxonomy = "references/silva_138_taxonomy.tsv",
    include.species = TRUE,
    compress = TRUE,
    n_euk = 0  # Exclude eukaryotes for prokaryote-focused studies
)

# Process UNITE database for fungi
buildUniteRef(
    fin = "sh_general_release_dynamic_25.07.2023.fasta",
    fout_fasta = "references/unite_formatted.fasta.gz", 
    fout_taxonomy = "references/unite_taxonomy.tsv",
    include.species = TRUE,
    compress = TRUE
)

# Verify outputs
silva_seqs <- Biostrings::readDNAStringSet("references/silva_138_formatted.fasta.gz")
silva_tax <- read.table("references/silva_138_taxonomy.tsv", sep="\t", header=TRUE)

cat("SILVA database contains", length(silva_seqs), "sequences\n")
cat("Taxonomy file contains", nrow(silva_tax), "records\n")
```

## Troubleshooting

### Common Issues

1.  **Memory Limitations**: For very large databases, consider increasing R memory limits or processing in batches

2.  **File Path Issues**: Ensure input files exist and output directories are writable

3.  **Format Compatibility**: Verify input files match expected SILVA/UNITE formats

### Error Messages

-   `"Input fasta file does not exist"`: Check file path and permissions

-   `"Duplicated sequence IDs detected"`: Input contains duplicate identifiers (automatically handled)

-   `"Failed to extract reference IDs"`: Input file format may not match expected structure

## Citation

If you use taxreftool in your research, please cite Callahan, B.J., et al. (2016) and:

```         
Suokas, M. (2025). taxreftool: Tools for Building Taxonomic Reference Databases. 
R package version 1.0.0. Retrieved from https://www.github.com/msuokas-bco/taxreftool
```

## License

This package is licensed under LGPL-3. See the LICENSE file for details.

## References

1.  Quast, C., et al. (2013). The SILVA ribosomal RNA gene database project. *Nucleic Acids Research*, 41(D1), D590-D596.

2.  Nilsson, R.H., et al. (2019). The UNITE database for molecular identification of fungi. *Nucleic Acids Research*, 47(D1), D259-D264.

3.  Callahan, B.J., et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*, 13(7), 581-583.
