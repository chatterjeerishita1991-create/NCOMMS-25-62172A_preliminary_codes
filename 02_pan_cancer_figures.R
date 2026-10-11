#!/usr/bin/env Rscript
# Figure 1A and 1B.

source("R/helpers.R")
check_packages()
cdr <- read_survival_data()
labels <- read_group_labels("pan_cancer")

## Figure 1A: TCGA-CDR overall survival with the checked-in R175/R273 labels.
clinical <- cdr |>
  dplyr::mutate(Mutation = labels$mutation[match(patient_barcode, labels$patient_barcode)]) |>
  dplyr::filter(!is.na(Mutation), !is.na(OS), !is.na(OS_time)) |>
  dplyr::mutate(Mutation = factor(Mutation, levels = c("R175", "R273")))

fit <- survival::survfit(survival::Surv(OS_time, OS) ~ Mutation, data = clinical)

# Keep this native to ggplot2/survival rather than relying on the ggpubr and
# survminer compatibility layer (which varies across Bioconductor containers).
surv_summary <- summary(fit, censored = TRUE)
curve_data <- data.frame(
  time = surv_summary$time,
  survival = surv_summary$surv,
  n_censor = surv_summary$n.censor,
  Mutation = sub("^Mutation=", "", as.character(surv_summary$strata))
)
curve_data <- dplyr::bind_rows(
  data.frame(time = 0, survival = 1, n_censor = 0, Mutation = levels(clinical$Mutation)),
  curve_data
) |>
  dplyr::mutate(Mutation = factor(Mutation, levels = levels(clinical$Mutation)))
logrank <- survival::survdiff(survival::Surv(OS_time, OS) ~ Mutation, data = clinical)
p_value <- stats::pchisq(logrank$chisq, df = length(logrank$n) - 1, lower.tail = FALSE)
time_breaks <- seq(0, 6000, by = 2000)
curve_plot <- ggplot2::ggplot(curve_data,
  ggplot2::aes(time, survival, color = Mutation)) +
  ggplot2::geom_step(linewidth = .7) +
  ggplot2::geom_point(data = dplyr::filter(curve_data, n_censor > 0), shape = 3, size = 1.8) +
  ggplot2::scale_color_manual(values = c(R175 = "#E64339", R273 = "#F79A36")) +
  ggplot2::scale_x_continuous(breaks = time_breaks, limits = range(time_breaks)) +
  ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25)) +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "Overall Survival by Mutation", x = "Time", y = "Survival probability", color = "Strata") +
  ggplot2::annotate("text", x = 900, y = .12, label = paste0("p = ", format.pval(p_value, digits = 2, eps = .001)))
risk_summary <- summary(fit, times = time_breaks, extend = TRUE)
risk_data <- expand.grid(time = time_breaks, Mutation = levels(clinical$Mutation)) |>
  dplyr::mutate(n_risk = as.vector(risk_summary$n.risk),
    Mutation = factor(Mutation, levels = levels(clinical$Mutation)))
risk_plot <- ggplot2::ggplot(risk_data, ggplot2::aes(time, Mutation, label = n_risk)) +
  ggplot2::geom_text(size = 3) +
  ggplot2::scale_x_continuous(breaks = time_breaks, limits = range(time_breaks)) +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.title.y = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size = 8)) +
  ggplot2::labs(title = "Number at risk", x = "Time")
grDevices::pdf(project_file("figures", "Figure_1A_survival.pdf"), width = 5, height = 6)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1,
  heights = grid::unit(c(4, 1.4), "null"))))
print(curve_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(risk_plot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
grDevices::dev.off()

## Figure 1B: VST CD274 in the exact five qualifying cancer types.
box_data <- lapply(c("UCEC", "STAD", "LGG", "COAD", "BRCA"), function(type) {
  path <- project_file("data", "cache", paste0("TCGA-", type, "_se.rds"))
  if (!file.exists(path)) stop("Run scripts/01_fetch_expression.R first: ", path)
  dds <- make_dds(readRDS(path))
  vst <- DESeq2::vst(dds, blind = TRUE)
  if (!"CD274" %in% rownames(vst)) stop("CD274 not present in ", type)
  tibble::tibble(Cancer = type, Mutation = as.character(vst$Mutation),
    Expression = SummarizedExperiment::assay(vst)["CD274", ]) |>
    dplyr::filter(Mutation %in% c("R175", "R273"))
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(Cancer = factor(Cancer, levels = c("UCEC", "STAD", "LGG", "COAD", "BRCA")),
    Mutation = factor(Mutation, levels = c("R175", "R273")))

star_label <- function(p) {
  if (is.na(p) || p > .05) "" else if (p <= .0001) "****" else if (p <= .001) "***" else if (p <= .01) "**" else "*"
}
comparison_labels <- box_data |>
  dplyr::group_by(Cancer) |>
  dplyr::summarise(
    p_value = tryCatch(stats::wilcox.test(Expression ~ Mutation)$p.value, error = function(e) NA_real_),
    y = max(Expression, na.rm = TRUE) + .07 * diff(range(Expression, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  dplyr::mutate(label = vapply(p_value, star_label, character(1)))
p_cd274 <- ggplot2::ggplot(box_data, ggplot2::aes(Cancer, Expression, color = Mutation)) +
  ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = .7), width = .6, outlier.shape = NA) +
  ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = .12, dodge.width = .7), alpha = .75, size = 1.2) +
  ggplot2::geom_text(data = dplyr::filter(comparison_labels, label != ""),
    ggplot2::aes(Cancer, y, label = label), inherit.aes = FALSE, color = "black", size = 5) +
  ggplot2::scale_color_manual(values = c(R175 = "#E64339", R273 = "#F79A36")) +
  ggplot2::theme_classic() + ggplot2::labs(title = "CD274", y = "VST Normalized Expression", x = "Cancer", color = "Mutation") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .2)))
save_plot(p_cd274, "Figure_1B_CD274.pdf", width = 6, height = 3)
utils::write.csv(box_data, project_file("results", "figure_1B_data.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), project_file("results", "session-pan-cancer.txt"))
