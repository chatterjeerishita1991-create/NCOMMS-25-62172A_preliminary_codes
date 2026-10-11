#!/usr/bin/env Rscript
# Figures 2G, 3E, S2G, S2H, S2I, and S3E: TCGA-OV R175H versus R273H.
# Plotting follows the CMCCTools functions used in the archived analysis.

source("R/helpers.R")
source("R/cmcc_plotting.R")
check_packages()

ov_path <- project_file("data", "cache", "TCGA-OV_se.rds")
if (!file.exists(ov_path)) stop("Run scripts/01_fetch_expression.R first: ", ov_path)
dds <- make_dds(readRDS(ov_path))
dds <- dds[, dds$Mutation %in% c("R175H", "R273H")]
dds$Mutation <- factor(dds$Mutation, levels = c("R273H", "R175H"))
DESeq2::design(dds) <- ~ Mutation
dds <- DESeq2::DESeq(dds)

# This is the same DESeq2 contrast and zero-p-value handling as the source Rmd.
marks <- DESeq2::results(dds, contrast = c("Mutation", "R175H", "R273H"))
rownames(marks) <- SummarizedExperiment::rowData(dds)$GENEID
marks$padj[marks$padj == 0] <- 1e-310
marks$pvalue[marks$pvalue == 0] <- 1e-310
de_export <- as.data.frame(marks) |>
  tibble::rownames_to_column("gene")
utils::write.csv(de_export, project_file("results", "ov_R175H_vs_R273H_DE.csv"), row.names = FALSE)

labels <- c(
  paste0("R175H n=", sum(dds$Mutation == "R175H")),
  paste0("R273H n=", sum(dds$Mutation == "R273H"))
)
main_volcano <- cmcc_gen_volcano(marks, labels, "Mutation R175H vs R273H",
  pCutoff = .05, subtitle = "All Samples")
save_plot(main_volcano[[1]], "Figure_S2G_volcano.pdf", 8, 8)

gsea_result <- cmcc_run_gsea(
  marks,
  project_file("data", "input", "gene_sets"),
  main_volcano[[2]], main_volcano[[3]], p_cutoff = .05
)
cmcc_gsea_csv(gsea_result$gsea, project_file("results", "ov_R175H_vs_R273H_GSEA.csv"))

save_cmcc_plot <- function(pathway, type, output) {
  key <- paste(pathway, type, sep = "::")
  plot <- gsea_result$plots[[key]]
  if (is.null(plot)) {
    stop("Required CMCCTools plot was not selected: ", pathway,
      ". Check the exact MSigDB 7.5.1 GMT files and cohort inputs.")
  }
  save_plot(plot, output, if (type == "enrichment") 7 else 8, if (type == "enrichment") 4 else 8)
}

save_cmcc_plot("Hm Epithelial Mesenchymal Transition", "enrichment", "Figure_2G_emt_enrichment.pdf")
save_cmcc_plot("Hm Tnfa Sig Via Nfkb", "enrichment", "Figure_3E_tnfa_nfkb_enrichment.pdf")
save_cmcc_plot("Hm Interferon Gamma Resp", "enrichment", "Figure_S3E_ifngamma_enrichment.pdf")
save_cmcc_plot("Hm Epithelial Mesenchymal Transition", "volcano", "Figure_S2H_emt_volcano.pdf")
save_cmcc_plot("Kg Oxidative Phosphorylation", "volcano", "Figure_S2I_kegg_oxphos_volcano.pdf")
writeLines(capture.output(sessionInfo()), project_file("results", "session-ovarian.txt"))
