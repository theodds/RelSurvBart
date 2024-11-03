additive_summary_plot_2 <- function (additive_summary, ribbonFill = "grey80", windsor = NA) {
  temp <- additive_summary$gamDf %>%
    mutate(term = case_when(
      term == "age" ~ "Age",
      term == "tpi" ~ "TPI",
      term == "wbc" ~ "WBC",
      TRUE           ~ term
    ))
  if (!is.na(windsor)) {
    if (!("quant" %in% colnames(temp))) {
      stop("Quantiles not supplied")
    }
    temp <- temp %>% filter(quant > windsor/2 & quant < 1 - windsor/2)
    glimpse(temp)
  }
  temp %>% distinct() %>%
    ggplot() + geom_hline(yintercept = 0) +
    geom_ribbon(aes(x_j, ymin = fx_j_lo, ymax = fx_j_hi), fill = ribbonFill, alpha = 0.5) +
    geom_line(aes(x_j, fx_j_mean), col = "firebrick3") +
    geom_rug(aes(x_j, fx_j_mean), sides = "b", alpha = 0.25) +
    facet_wrap(~term, scale = "free") +
    labs(x = ("term"), y = ("Partial effect"))
}
