## =====================================================================
## simulate_data.R  —  Synthetic inputs for the micro-simulation showcase
## ---------------------------------------------------------------------
## Everything here is 100% synthetic and generated on the fly: no real
## data of any kind. It produces the same *shape* of inputs a health
## micro-simulation needs, so the engines have something realistic to run.
## =====================================================================

## The toy model: a synthetic population of individuals who, year by year,
## can develop chronic diseases and eventually die. Each individual carries
## a per-disease "exposure" relative risk (a stand-in for things like air
## pollution or physical-activity effects). Baseline incidence/mortality
## rates depend on age and sex. A simple comorbidity rule makes some
## diseases raise the risk of others.

DISEASES <- c("asthma", "diabetes", "heart_disease", "depression",
              "copd", "stroke", "all_cause_mortality")
CHRONIC  <- setdiff(DISEASES, "all_cause_mortality")

## --- baseline per-year rates by age band and sex (made up, plausible) ----
## returned as an array [age 0..100, sex 1..2, disease] for O(1) lookup
build_rates <- function(seed = 1) {
  set.seed(seed)
  ages <- 0:100
  arr <- array(0, dim = c(length(ages), 2, length(DISEASES)),
               dimnames = list(age = ages, sex = c("male", "female"), disease = DISEASES))
  for (d in DISEASES) {
    base <- switch(d,
      asthma           = 0.004, diabetes      = 0.006, heart_disease = 0.005,
      depression       = 0.008, copd          = 0.003, stroke        = 0.004,
      all_cause_mortality = 0.006)
    for (s in 1:2) {
      ## rates rise with age; mortality rises steeply
      slope <- if (d == "all_cause_mortality") 0.0016 else 0.00035
      male_factor <- if (s == 1) 1.15 else 0.9
      arr[, s, d] <- pmax(0, base * (1 + slope * ages * 100) * male_factor)
    }
  }
  arr
}

## --- comorbidity relative risks (which disease raises which) --------------
## one row per (risk_factor disease -> outcome disease), multiplicative
build_comorbidity <- function() {
  data.frame(
    risk_factor = c("diabetes", "copd",  "depression", "heart_disease", "diabetes"),
    outcome     = c("heart_disease", "heart_disease", "stroke", "stroke", "stroke"),
    rr          = c(1.8, 1.5, 1.4, 2.0, 1.6),
    stringsAsFactors = FALSE
  )
}

## --- the synthetic population --------------------------------------------
simulate_population <- function(n = 20000, seed = 42) {
  set.seed(seed)
  age <- pmin(95L, pmax(18L, as.integer(round(rnorm(n, 45, 16)))))
  sex <- sample(1:2, n, replace = TRUE)            # 1 male, 2 female
  pop <- data.frame(id = seq_len(n), age = age, sex = sex)
  ## per-individual exposure RR for each disease (log-normal, mean ~1)
  for (d in DISEASES) pop[[d]] <- exp(rnorm(n, 0, 0.15))
  ## baseline prevalence: a few % already have each chronic disease at start
  init <- matrix(FALSE, n, length(CHRONIC), dimnames = list(NULL, CHRONIC))
  for (d in CHRONIC) init[, d] <- runif(n) < 0.03
  list(pop = pop, init = init)
}

## --- convenience: assemble the full input bundle --------------------------
make_inputs <- function(n = 20000, n_cycles = 20, seed = 42) {
  dat <- simulate_population(n, seed)
  list(pop        = dat$pop,
       init       = dat$init,
       rates      = build_rates(seed),
       comorb     = build_comorbidity(),
       diseases   = DISEASES,
       chronic    = CHRONIC,
       n_cycles   = n_cycles)
}
