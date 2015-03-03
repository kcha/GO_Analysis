# Jan 23, 2013
#
# Perform GO enrichment

#' @export
load_genes <- function(file, db) {
    # Input is single column of gene symbols
    write(paste("// Testing", file), stderr())
    zz <- read.table(file, header=F, sep="\t")
    zzz <- db[db$alias_symbol %in% zz[,1],]
    return(zzz)
}

#' @export
load_universe <- function(file, db) {
    # Universe of genes
    #   Input is single column of gene symbols
    symbols <- read.table(file, header=F, sep="\t")
    symbolsz <- db[db$alias_symbol %in% symbols[,1],]
    return(symbolsz)
}

#' @export
load_alias2eg <- function() {
    return(toTable(org.Hs.egALIAS2EG))
}

#' @export
go_test <- function(genes, universeGenes, ont = "BP") {
    params <- new("GOHyperGParams",
                  geneIds = genes,
                  universeGeneIds = universeGenes,
                  annotation = "org.Hs.eg.db",
                  ontology = ont,
                  pvalueCutoff = 0.05,
                  testDirection = "over",
                  conditional = TRUE)
    
    hgOver <- hyperGTest(params)
    
    return(hgOver)
}

#' @export
get_sig_terms <- function(hgOver, filtered) {
    s <- summary(hgOver, pvalue=1)
    
    pvals.adj <- p.adjust(s[,2], method = "BH")
    
    s2 <- cbind(s, pvals.adj)
    
    colnames(s2)[ncol(s2)] <- "AdjustedPValue"
    
    # sort results by adjusted p-value
    s2 <- s2[order(s2[,ncol(s2)]),]
    
    return(s2)
}

#' @export
go_enrichment_2_table <- function(result) {
    # Convert significant hits to a table
    # Combine all three ontologies
    Ontologies <- c("BP", "MF", "CC")
    SigTermsTab <- vector("list", length=length(Ontologies))
    
    for (ont in c("BP", "MF", "CC")) {
        # Check that there are significant GO terms
        if (nrow(result[[ont]]$summary) > 0) {
            SigTermsTab[[ont]] <- cbind(result[[ont]]$summary, Ontology=ont)
            colnames(SigTermsTab[[ont]])[1] <- "GOID"
        }
    }

    ALL <- do.call("rbind", SigTermsTab)
    return(ALL)
}

#' @export
write_results <- function(g, dir, inputfile) {
    tab <- go_enrichment_2_table(g)
    #output <- paste("go_enrichment", "tab", sep=".")
    outfile <- sub("\\.[^.]+$", ".go_enrichment.tab", inputfile)
    write(paste("// Writing to", outfile), stderr())
    write.table(tab, file=outfile, quote=F, sep="\t", row.names=FALSE)
    
    # Write results in Enrichment Map format
    # The Generic Enrichment Results file needs:
    #   gene-set ID (must match the gene-set ID in the GMT file),
    #   gene-set name or description,
    #   p-value,
    #   FDR correction value
    #   Phenotype: +1 or -1, to identify enrichment in up- and down-regulation, or, more in general, in either of the two phenotypes being compared in the two-class analysis
    #       +1 maps to red
    #       -1 maps to blue
    #if (!is.null(tab)) {
        #tab <- tab[, c(1,7,2,8)]
        #tab <- cbind(tab, Phenotype=1)
        #output <- paste("EM_EnrTable", clustno, "Generic", "txt", sep=".")
        #write.table(tab, file=file.path(dir, output), quote=F, sep="\t", 
                    #row.names=FALSE, col.names=TRUE)
    #}
}

#' @export
run_go_test <- function(file, universe, db = NULL) {
    if (is.null(db))
        db <- load_alias2eg()
    
    G <- load_genes(file, db)
    U <- load_universe(universe, db)
    
    hg.bp <- go_test(G[,1], U[,1], ont = "BP")
    hg.mf <- go_test(G[,1], U[,1], ont = "MF")
    hg.cc <- go_test(G[,1], U[,1], ont = "CC")
    
    pv.bp <- get_sig_terms(hg.bp)
    pv.mf <- get_sig_terms(hg.mf)
    pv.cc <- get_sig_terms(hg.cc)
    
    BP <- list(hyperg=hg.bp, summary=pv.bp)
    MF <- list(hyperg=hg.mf, summary=pv.mf)
    CC <- list(hyperg=hg.cc, summary=pv.cc)
    
    return(list(BP=BP, MF=MF, CC=CC))
}
