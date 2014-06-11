#!/usr/bin/env Rscript
# Feb 6, 2013
#
# Perform GO enrichment on all clusters

print_help <- function() {
  text <- "**** GO Enrichment Analysis in R ****
Script for performing GO enrichment analysis in R

Usage: ./run_go_test.R universe_genes.txt target_genes.txt [target_genes_2.txt ...]

Arguments:
  1) Universe genes - a list of universe genes to be tested against during GO 
    analysis. Requires gene symbols in column 1.
  2) Target genes (one or more) - a list of genes to be tested for GO enrichment. 
    Requires gene symbols in column 1.

Options:
  --help        This help message

Output:
  A tab-delimited file of enriched GO terms with statistics for each target file.

Test run:
  ./run_go_test.R test_data/GENE_UNIVERSE.tab test_data/sample_genes.txt
"
  writeLines(text, stderr())
}


args <- commandArgs(TRUE)

if (length(args) < 1) {
  print_help()
  stop("Missing arguments")
}
if (args[1] %in% c("-h", "--help", "-help")) {
  print_help()
  stop("Terminating")
}

library(stringr)
source("go_test_func.r")

universe <- args[1]
files <- args[2:length(args)]

###############################################################################

write("// Running enrichment...", stderr())
db <- load_alias2eg()
#go_enrich <- lapply(file.path(dir, files), run_go_test, db = db)
go_enrich <- lapply(files, run_go_test, universe = universe, db = db)
names(go_enrich) <- files

####
#savefile <- paste(id, getDate(TRUE), "cluster.go_enrichment", "RData", sep=".")
#save(files, db, go_enrich, file=file.path(dir, savefile), compress=T)
#print(paste("Saving results to", savefile))

####
# Print results
z <- lapply(1:length(go_enrich), function(x) write_results(go_enrich[[x]], dir, files[x]))

####
write("// Done!", stderr())
