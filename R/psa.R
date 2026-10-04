## =====================================================================
## psa.R  —  Uncertainty layer: Latin Hypercube Sampling + PSA
## ---------------------------------------------------------------------
## Probabilistic Sensitivity Analysis (PSA) runs the model many times,
## each time drawing the uncertain parameters from their distributions,
## to turn a single point estimate into a distribution with credible
## intervals. Latin Hypercube Sampling (LHS) picks those draws in a
## stratified way, so we cover the parameter space well with FAR fewer
## runs than plain Monte Carlo.
##
## Here the uncertain parameters are the comorbidity relative risks:
## each has a central value and a (synthetic) confidence interval.
## =====================================================================

## --- Latin Hypercube Sampling: n_draws points in [0,1]^k -----------------
lhs_unit <- function(n_draws, k, seed = 1) {
  set.seed(seed)
  u <- matrix(0, n_draws, k)
  for (j in seq_len(k)) {
    ## one stratified point per interval [i/n, (i+1)/n), then shuffle
    strata <- (sample(n_draws) - 1 + runif(n_draws)) / n_draws
    u[, j] <- strata
  }
  u
}

## --- turn unit draws into log-normal RR draws from a 95% CI ---------------
## given central rr and ci = c(low, high), treat as log-normal
draw_rr_from_ci <- function(u_col, central, ci) {
  logmean <- log(central)
  logsd   <- (log(ci[2]) - log(ci[1])) / (2 * 1.96)
  exp(qnorm(u_col, logmean, logsd))
}

## --- build the PSA draw table for the comorbidity RRs ---------------------
## comorb_ci: data.frame with risk_factor, outcome, rr, ci_low, ci_high
build_psa_draws <- function(comorb_ci, n_draws = 200, seed = 1) {
  k <- nrow(comorb_ci)
  u <- lhs_unit(n_draws, k, seed)
  draws <- matrix(0, n_draws, k)
  for (j in seq_len(k)) {
    draws[, j] <- draw_rr_from_ci(u[, j], comorb_ci$rr[j],
                                  c(comorb_ci$ci_low[j], comorb_ci$ci_high[j]))
  }
  colnames(draws) <- paste(comorb_ci$risk_factor, comorb_ci$outcome, sep = "->")
  draws
}

## --- run PSA: for each draw, overwrite RRs, run the fast engine, reduce ----
## returns an array [n_draws, n_cycles, n_chronic] of prevalence — but we
## keep only the per-cycle prevalence (already reduced inside the engine),
## never the full individual state, so memory stays flat.
run_psa <- function(inp, comorb_ci, n_draws = 200, seed = 1,
                    engine = run_engine_fast, progress = TRUE) {
  draws <- build_psa_draws(comorb_ci, n_draws, seed)
  chronic <- inp$chronic
  out <- array(NA_real_, dim = c(n_draws, inp$n_cycles, length(chronic)),
               dimnames = list(draw = NULL, cycle = NULL, disease = chronic))
  for (i in seq_len(n_draws)) {
    inp_i <- inp
    inp_i$comorb <- data.frame(risk_factor = comorb_ci$risk_factor,
                               outcome     = comorb_ci$outcome,
                               rr          = draws[i, ],
                               stringsAsFactors = FALSE)
    res <- engine(inp_i, seed = 1000 + i)   # reduced result only
    out[i, , ] <- res$prevalence
    if (progress && i %% max(1, n_draws %/% 10) == 0)
      message(sprintf("  PSA draw %d/%d", i, n_draws))
  }
  out
}

## --- summarise PSA output into median + 95% credible interval ------------
summarise_psa <- function(psa_out) {
  chronic <- dimnames(psa_out)$disease
  n_cyc   <- dim(psa_out)[2]
  rows <- list()
  for (j in seq_along(chronic)) {
    for (cyc in seq_len(n_cyc)) {
      v <- psa_out[, cyc, j]
      rows[[length(rows) + 1]] <- data.frame(
        disease = chronic[j], cycle = cyc,
        median = median(v),
        lo95 = unname(quantile(v, 0.025)),
        hi95 = unname(quantile(v, 0.975)))
    }
  }
  do.call(rbind, rows)
}

## --- the synthetic comorbidity CIs (central value + made-up 95% CI) -------
default_comorb_ci <- function() {
  cb <- build_comorbidity()
  cb$ci_low  <- cb$rr / 1.25
  cb$ci_high <- cb$rr * 1.25
  cb
}
