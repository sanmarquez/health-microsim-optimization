## =====================================================================
## engine_fast.R  —  The "after" engine (vectorized, ~100x faster)
## ---------------------------------------------------------------------
## Same model, same transition rule, same result distribution — but the
## state is a LOGICAL MATRIX H (individuals x diseases) instead of text,
## and every operation works on whole columns at once. No per-person
## sapply, no string parsing. This is the core optimization the showcase
## is about.
##
## Transition rule (identical to the naive engine):
##     p = 1 - exp(-(rate * exposure_RR * comorbidity_RR))
## =====================================================================

run_engine_fast <- function(inp, seed = 1) {
  set.seed(seed)
  pop      <- inp$pop
  n        <- nrow(pop)
  diseases <- inp$diseases
  chronic  <- inp$chronic
  D        <- length(chronic)
  rates    <- inp$rates
  comorb   <- inp$comorb
  sex      <- pop$sex

  ## state: logical matrix, one column per chronic disease (+ alive flag)
  H <- inp$init[, chronic, drop = FALSE]        # n x D logical
  storage.mode(H) <- "logical"
  alive <- rep(TRUE, n)
  col_of <- setNames(seq_len(D), chronic)

  ## precompute comorbidity as (outcome -> list of {rr, risk_col})
  comorb_by_out <- split(comorb, comorb$outcome)

  prevalence <- matrix(0, inp$n_cycles, D, dimnames = list(NULL, chronic))
  deaths <- integer(inp$n_cycles)

  for (cyc in seq_len(inp$n_cycles)) {
    age <- pmin(100L, pop$age + cyc - 1L)
    age_idx <- age + 1L

    for (d in diseases) {
      ## base rate: ONE vectorized array lookup for everyone
      base <- rates[cbind(age_idx, sex, match(d, diseases))]
      expo <- pop[[d]]

      ## comorbidity RR: start at 1, multiply in each relevant disease column
      cf <- rep(1, n)
      cb <- comorb_by_out[[d]]
      if (!is.null(cb)) {
        for (r in seq_len(nrow(cb))) {
          rc <- col_of[[cb$risk_factor[r]]]
          if (!is.null(rc)) cf <- cf * ifelse(H[, rc], cb$rr[r], 1)
        }
      }

      p <- 1 - exp(-(base * expo * cf))
      u <- runif(n)
      fires <- alive & (u < p) & (age < 100)

      if (d == "all_cause_mortality") {
        alive[fires] <- FALSE
      } else {
        H[fires, col_of[[d]]] <- TRUE          # whole-column assignment
      }
    }

    ## tally: column sums over the living — no string parsing
    live <- alive
    prevalence[cyc, ] <- colSums(H[live, , drop = FALSE])
    deaths[cyc] <- sum(!alive)
  }

  list(prevalence = prevalence, deaths = deaths)
}
