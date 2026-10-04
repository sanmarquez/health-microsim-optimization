## =====================================================================
## engine_naive.R  —  The "before" engine (readable but slow)
## ---------------------------------------------------------------------
## This mirrors how a micro-simulation is often written first: each
## individual's health state is a STRING ("asthma depression" / "healthy"
## / "dead"), and every cycle we loop person-by-person, splitting and
## pasting text and calling sapply() to compute comorbidity risk.
##
## It is correct and easy to read — and it is the bottleneck the fast
## engine removes. Transition rule (shared by both engines):
##     p = 1 - exp(-(rate * exposure_RR * comorbidity_RR))
## =====================================================================

## comorbidity multiplier for one individual, given their disease string
.naive_comorb_rr <- function(dis_string, outcome, comorb) {
  present <- strsplit(dis_string, " ", fixed = TRUE)[[1]]
  rows <- comorb[comorb$outcome == outcome & comorb$risk_factor %in% present, ]
  if (nrow(rows) == 0) return(1)
  prod(rows$rr)
}

run_engine_naive <- function(inp, seed = 1) {
  set.seed(seed)
  pop      <- inp$pop
  n        <- nrow(pop)
  diseases <- inp$diseases
  chronic  <- inp$chronic
  rates    <- inp$rates
  comorb   <- inp$comorb

  ## state as a character vector; seed initial prevalence into the strings
  state <- rep("healthy", n)
  for (i in seq_len(n)) {
    have <- chronic[inp$init[i, ]]
    if (length(have)) state[i] <- paste(have, collapse = " ")
  }

  prevalence <- matrix(0, inp$n_cycles, length(chronic),
                       dimnames = list(NULL, chronic))
  deaths <- integer(inp$n_cycles)

  for (cyc in seq_len(inp$n_cycles)) {
    age <- pmin(100L, pop$age + cyc - 1L)

    for (d in diseases) {
      ## base rate for each individual (string-key-free, but per-person)
      base <- vapply(seq_len(n), function(i) rates[age[i] + 1L, pop$sex[i], d], numeric(1))
      expo <- pop[[d]]

      ## comorbidity RR, person by person, by parsing the state string
      cf <- vapply(seq_len(n), function(i) {
        if (state[i] == "dead") return(1)
        .naive_comorb_rr(state[i], d, comorb)
      }, numeric(1))

      p <- 1 - exp(-(base * expo * cf))
      u <- runif(n)

      alive <- state != "dead"
      fires <- alive & (u < p) & (age < 100)

      if (d == "all_cause_mortality") {
        state[fires] <- "dead"
      } else {
        for (i in which(fires)) {
          if (!grepl(d, state[i], fixed = TRUE)) {
            state[i] <- if (state[i] == "healthy") d else paste(state[i], d)
          }
        }
      }
    }

    ## tally this cycle by parsing strings
    for (j in seq_along(chronic)) {
      dd <- chronic[j]
      prevalence[cyc, j] <- sum(vapply(state, function(s)
        s != "dead" && grepl(dd, s, fixed = TRUE), logical(1)))
    }
    deaths[cyc] <- sum(state == "dead")
  }

  list(prevalence = prevalence, deaths = deaths)
}
