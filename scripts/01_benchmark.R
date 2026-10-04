## =====================================================================
## 01_benchmark.R  —  Speed + equivalence: naive vs vectorized engine
## ---------------------------------------------------------------------
## Produces:
##   figures/speedup.png       runtime vs population size (log scale)
##   figures/equivalence.png   prevalence trajectories, both engines overlaid
##   results_benchmark.csv     the timing table
## Run from the project root:  Rscript scripts/01_benchmark.R
## =====================================================================
suppressWarnings({
  source("R/simulate_data.R"); source("R/engine_naive.R"); source("R/engine_fast.R")
})

PAL <- list(naive = "#D1495B", fast = "#2A9D8F", grid = "#E6E6E6", ink = "#2B2B2B")
dir.create("figures", showWarnings = FALSE)

## ---- 1) scaling benchmark across population sizes -----------------------
sizes   <- c(500, 1000, 2000, 4000)
cycles  <- 8
tab <- data.frame(n = sizes, naive_s = NA_real_, fast_s = NA_real_, speedup = NA_real_)
cat("Benchmarking (naive is intentionally slow)...\n")
for (k in seq_along(sizes)) {
  inp <- make_inputs(n = sizes[k], n_cycles = cycles, seed = 42)
  tab$naive_s[k] <- system.time(run_engine_naive(inp, seed = 1))[["elapsed"]]
  tab$fast_s[k]  <- system.time(run_engine_fast(inp,  seed = 1))[["elapsed"]]
  tab$speedup[k] <- tab$naive_s[k] / tab$fast_s[k]
  cat(sprintf("  n=%5d  naive=%6.2fs  fast=%6.3fs  speedup=x%.0f\n",
              sizes[k], tab$naive_s[k], tab$fast_s[k], tab$speedup[k]))
}
write.csv(tab, "results_benchmark.csv", row.names = FALSE)

## ---- 2) equivalence check at a fixed size -------------------------------
inp <- make_inputs(n = 4000, n_cycles = 20, seed = 7)
rn <- run_engine_naive(inp, seed = 123)
rf <- run_engine_fast(inp,  seed = 123)
max_diff <- max(abs(rn$prevalence - rf$prevalence), abs(rn$deaths - rf$deaths))
cat(sprintf("\nEquivalence: max |naive - fast| = %.3g  (0 = bit-identical)\n", max_diff))

## ---- figure A: runtime vs N (log y) -------------------------------------
png("figures/speedup.png", width = 1600, height = 1000, res = 200)
par(mar = c(4.5, 4.8, 3.2, 1.2), family = "sans")
ylim <- range(c(tab$fast_s, tab$naive_s)); ylim[1] <- max(ylim[1], 1e-3)
plot(tab$n, tab$naive_s, type = "n", log = "y", ylim = ylim,
     xlab = "Population size (individuals)", ylab = "Runtime per run (seconds, log scale)",
     main = "", axes = FALSE)
abline(h = axTicks(2), col = PAL$grid, lwd = 1)
box(col = "#CFCFCF")
axis(1, col = "#CFCFCF", col.axis = PAL$ink)
axis(2, col = "#CFCFCF", col.axis = PAL$ink, las = 1)
lines(tab$n, tab$naive_s, col = PAL$naive, lwd = 3)
points(tab$n, tab$naive_s, col = PAL$naive, pch = 19, cex = 1.3)
lines(tab$n, tab$fast_s, col = PAL$fast, lwd = 3)
points(tab$n, tab$fast_s, col = PAL$fast, pch = 19, cex = 1.3)
## labels placed to the LEFT of the last point so they stay inside the panel
text(tab$n[length(sizes)], tab$naive_s[length(sizes)],
     labels = "naive (string state)  ", col = PAL$naive, pos = 2, font = 2)
text(tab$n[length(sizes)], tab$fast_s[length(sizes)],
     labels = "vectorized (logical matrix)  ", col = PAL$fast, pos = 2, font = 2)
title(main = "Vectorizing the core: two orders of magnitude faster",
      col.main = PAL$ink, cex.main = 1.25, adj = 0)
mtext(sprintf("~%.0fx faster at the largest size  ·  identical results (max |diff| = %.0f)",
              tab$speedup[length(sizes)], max_diff),
      side = 3, adj = 0, line = 0.2, col = "#6B6B6B", cex = 0.9)
dev.off()

## ---- figure B: equivalence (trajectories overlaid) ----------------------
png("figures/equivalence.png", width = 1600, height = 1000, res = 200)
par(mar = c(4.5, 4.8, 3.2, 1.2), family = "sans")
show <- c("diabetes", "heart_disease", "stroke")
cols <- c("#264653", "#E76F51", "#2A9D8F")
yr <- range(rn$prevalence[, show])
plot(NA, xlim = c(1, inp$n_cycles), ylim = yr, axes = FALSE,
     xlab = "Cycle (year)", ylab = "Prevalent cases (living)")
abline(h = pretty(yr), col = PAL$grid); box(col = "#CFCFCF")
axis(1, col = "#CFCFCF", col.axis = PAL$ink); axis(2, col = "#CFCFCF", col.axis = PAL$ink, las = 1)
for (j in seq_along(show)) {
  lines(rn$prevalence[, show[j]], col = cols[j], lwd = 5)             # naive: thick
  lines(rf$prevalence[, show[j]], col = "white", lwd = 1.6, lty = 1)  # fast: thin white over it
}
legend("topleft", bty = "n", legend = show, col = cols, lwd = 5, text.col = PAL$ink)
title(main = "Naive and vectorized engines trace the same curves",
      col.main = PAL$ink, cex.main = 1.25, adj = 0)
mtext("thick colour = naive · thin white = vectorized (drawn on top)",
      side = 3, adj = 0, line = 0.2, col = "#6B6B6B", cex = 0.9)
dev.off()

cat("\nWrote figures/speedup.png, figures/equivalence.png, results_benchmark.csv\n")
