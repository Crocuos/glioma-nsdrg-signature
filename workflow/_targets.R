suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(survival)
  library(ggplot2); library(timeROC); library(survminer); library(rms)
})


load_cohort <- function(expr_file, clin_file){
  expr <- readRDS(expr_file)           # genes x samples (log-normalized)
  clin <- fread(clin_file) |> as.data.frame()
  rownames(clin) <- clin$sample
  common <- intersect(colnames(expr), rownames(clin))
  list(expr = expr[, common, drop = FALSE], clin = clin[common, , drop = FALSE])
}


score_signature <- function(expr, coef_vec){
  g <- intersect(names(coef_vec), rownames(expr))
  if(length(g) < length(coef_vec))
    warning("Missing genes: ", paste(setdiff(names(coef_vec), g), collapse=", "))
  risk <- colSums(t(expr[g, , drop = FALSE]) * coef_vec[g])
  data.frame(sample = names(risk), risk = as.numeric(risk), row.names = names(risk))
}


km_plot <- function(clin, risk, cutpoint, label){
  df <- data.frame(clin[rownames(clin), c("os_time","os_event")],
                   risk = risk[rownames(clin), "risk"])
  df$group <- ifelse(df$risk >= cutpoint, "High", "Low")
  fit <- survfit(Surv(os_time, os_event) ~ group, data = df)
  p <- ggsurvplot(fit, data = df, risk.table = FALSE, pval = TRUE,
                  legend.title = label, legend.labs = c("High","Low"),
                  palette = c("#D55E00","#0072B2"))
  out <- file.path("results", paste0("KM_", label, ".pdf"))
  dir.create("results", showWarnings = FALSE)
  ggsave(out, p$plot, width = 4.5, height = 4)
  out
}


roc_plot <- function(clin, risk, label, times = c(365, 1095, 1825)){
  df <- data.frame(time = clin$os_time, status = clin$os_event,
                   risk = risk[rownames(clin), "risk"])
  roc <- timeROC(T = df$time, delta = df$status, marker = df$risk,
                 cause = 1, times = times, iid = TRUE)
  dd <- data.frame(t = roc$times, AUC = roc$AUC)
  p <- ggplot(dd, aes(t/365, AUC)) + geom_line() + geom_point() +
       labs(x="Years", y="AUC", title=paste0("timeROC - ", label)) +
       ylim(0.5, 1)
  out <- file.path("results", paste0("timeROC_", label, ".pdf"))
  dir.create("results", showWarnings = FALSE); ggsave(out, p, width = 4.5, height = 4)
  out
}


calib_plot <- function(clin, risk, label, years = c(1,3,5)){
  df <- data.frame(os_time = clin$os_time, os_event = clin$os_event,
                   risk = risk[rownames(clin), "risk"])
  dd <- datadist(df); options(datadist = "dd")
  fit <- cph(Surv(os_time, os_event) ~ risk, data = df, x = TRUE, y = TRUE, surv = TRUE)
  S <- Survival(fit)
  pred <- sapply(years, function(y) 1 - S(y*365, lp = df$risk))
  colnames(pred) <- paste0(years, "y")

  library(tidyr)
  long <- cbind(df, pred) |> pivot_longer(starts_with(c("1y","3y","5y")),
                                          names_to="h", values_to="phat")
  out <- file.path("results", paste0("calibration_", label, ".pdf"))
  p <- ggplot(long, aes(phat, ..scaled..)) + geom_histogram(bins=20) +
       facet_wrap(~h, nrow=1) + labs(x="Predicted risk", y="Scaled density",
       title=paste0("Calibration (qualitative) - ", label))
  dir.create("results", showWarnings = FALSE); ggsave(out, p, width = 7, height = 3)
  out
}

save_coeffs <- function(coef_vec, path){
  dir.create(dirname(path), showWarnings = FALSE)
  write.csv(data.frame(Gene = names(coef_vec), Coef = coef_vec), path, row.names = FALSE)
  path
}
write_session <- function(path){
  dir.create(dirname(path), showWarnings = FALSE)
  writeLines(c(capture.output(sessionInfo())), path); path
}

