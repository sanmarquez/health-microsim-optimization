## =====================================================================
## 02_psa.R  —  Uncertainty: Latin Hypercube Sampling + PSA intervals
## ---------------------------------------------------------------------
## Runs the fast engine many times, each with the comorbidity relative
## risks drawn (via LHS) from their confidence intervals, and turns the
## point estimate into a median + 95% credible band.
## Produces:
##   figures/psa_intervals.png   prevalence with uncertainty band
##   figures/lhs_coverage.png    LHS vs plain Monte Carlo coverage
##   results_psa.csv             median / lo95 / hi95 per cycle & disease
## Run from the project root:  Rscript scripts/02_psa.R
## =====================================================================
suppressWarnings({
  source("R/simulate_data.R"); source("R/engine_fast.R"); source("R/psa.R")
})
PAL <- list(band = "#9AD1C9", line = "#2A9D8F", ink = "#2B2B2B", grid = "#E6E6E6",
            lhs = "#2A9D8F", mc = "#D1495B")
dir.create("figures", showWarnings = FALSE)

N_DRAWS <- 200
inp     <- make_inputs(n = 5000, n_cycles = 20, seed = 42)
cb_ci   <- default_comorb_ci()

cat(sprintf("Running PSA: %d LHS draws x fast engine...\n", N_DRAWS))
psa_out <- run_psa(inp, cb_ci, n_draws = N_DRAWS, seed = 1, progress = TRUE)
summ    <- summarise_psa(psa_out)
write.csv(summ, "results_psa.csv", row.names = FALSE)

## ---- figure A: a prevalence trajectory with its 95% credible band -------
target <- "heart_disease"
s <- summ[summ$disease == target, ]
png("figures/psa_intervals.png", width = 1600, height = 1000, res = 200)
par(mar = c(4.5, 4.8, 3.2, 1.2), family = "sans")
plot(NA, xlim = range(s$cycle), ylim = range(c(s$lo95, s$hi95)), axes = FALSE,
     xlab = "Cycle (year)", ylab = sprintf("Prevalent cases of %s", target))
abline(h = pretty(range(c(s$lo95, s$hi95))), col = PAL$grid); box(col = "#CFCFCF")
axis(1, col = "#CFCFCF", col.axis = PAL$ink); axis(2, col = "#CFCFCF", col.axis = PAL$ink, las = 1)
polygon(c(s$cycle, rev(s$cycle)), c(s$lo95, rev(s$hi95)),
        col = PAL$band, border = NA)
lines(s$cycle, s$median, col = PAL$line, lwd = 3)
legend("topleft", bty = "n",
       legend = c("median", "95% credible interval"),
       col = c(PAL$line, PAL$band), lwd = c(3, 10), text.col = PAL$ink)
title(main = "Parameter uncertainty propagated to the outcome",
      col.main = PAL$ink, cex.main = 1.25, adj = 0)
mtext(sprintf("%d LHS draws over the comorbidity relative risks", N_DRAWS),
      side = 3, adj = 0, line = 0.2, col = "#6B6B6B", cex = 0.9)
dev.off()

## ---- figure B: why LHS? coverage of the unit interval vs plain MC --------
set.seed(1)
nshow <- 20
lhs_pts <- lhs_unit(nshow, 1, seed = 3)[, 1]
mc_pts  <- runif(nshow)
png("figures/lhs_coverage.png", width = 1600, height = 700, res = 200)
par(mar = c(4.2, 6.5, 3.0, 1.2), family = "sans")
plot(NA, xlim = c(0, 1), ylim = c(0.5, 2.5), axes = FALSE, xlab = "Parameter quantile (0-1)", ylab = "")
abline(v = seq(0, 1, by = 1/nshow), col = PAL$grid)
points(lhs_pts, rep(2, nshow), pch = 19, col = PAL$lhs, cex = 1.3)
points(mc_pts,  rep(1, nshow), pch = 19, col = PAL$mc,  cex = 1.3)
axis(1, col = "#CFCFCF", col.axis = PAL$ink)
mtext(c("plain\nMonte Carlo", "Latin\nHypercube"), side = 2, at = c(1, 2),
      las = 1, line = 1, col = c(PAL$mc, PAL$lhs), font = 2)
title(main = sprintf("LHS covers the space evenly with the same %d draws", nshow),
      col.main = PAL$ink, cex.main = 1.2, adj = 0)
mtext("each vertical band should hold one point; MC clumps and leaves gaps",
      side = 3, adj = 0, line = 0.2, col = "#6B6B6B", cex = 0.85)
dev.off()

## ---- console summary ----------------------------------------------------
final <- summ[summ$cycle == inp$n_cycles, ]
final$width_rel <- round((final$hi95 - final$lo95) / pmax(final$median, 1), 3)
cat("\nFinal-cycle PSA summary (median [95% CI]):\n")
for (i in seq_len(nrow(final)))
  cat(sprintf("  %-14s %6.0f  [%6.0f, %6.0f]  rel.width %.2f\n",
              final$disease[i], final$median[i], final$lo95[i], final$hi95[i], final$width_rel[i]))
cat("\nWrote figures/psa_intervals.png, figures/lhs_coverage.png, results_psa.csv\n")
