#!/usr/bin/env Rscript
# Downloads the smallest STAR-count cohort needed by the requested figures.

source("R/helpers.R")
check_packages()
dir.create(project_file("data", "cache"), recursive = TRUE, showWarnings = FALSE)

pan_labels <- read_group_labels("pan_cancer")
ov_labels <- read_group_labels("ovarian")

# These are exactly the five disease types shown in Figure 1B. OV is used only
# for the R175H-versus-R273H DE/GSEA figures.
for (type in c("BRCA", "COAD", "LGG", "STAD", "UCEC")) {
  output <- project_file("data", "cache", paste0("TCGA-", type, "_se.rds"))
  if (!file.exists(output)) {
    message("Downloading TCGA-", type)
    se <- prepare_star_counts(paste0("TCGA-", type), project_file("data", "cache"))
    se <- apply_group_labels(se, pan_labels)
    saveRDS(se, output)
  }
}

ov_output <- project_file("data", "cache", "TCGA-OV_se.rds")
if (!file.exists(ov_output)) {
  message("Downloading TCGA-OV")
  ov <- prepare_star_counts("TCGA-OV", project_file("data", "cache"))
  ov <- apply_group_labels(ov, ov_labels)
  saveRDS(ov, ov_output)
}

writeLines(capture.output(sessionInfo()), project_file("results", "session-fetch.txt"))
