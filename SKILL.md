---
name: run-pipeline
description: Run, reproduce, or extend this micro-simulation showcase. Use when the user asks to run the benchmark or PSA, regenerate the figures, render the Quarto report, change the population size / number of PSA draws, or add a disease or comorbidity.
---

# Running the micro-simulation showcase

This repo is a self-contained R demo. Base R only; no packages required for the
two scripts. Quarto is only needed for the HTML report.

## Reproduce everything
```bash
Rscript scripts/01_benchmark.R   # speed + equivalence -> figures/, results_benchmark.csv
Rscript scripts/02_psa.R         # LHS/PSA intervals    -> figures/, results_psa.csv
quarto render report.qmd         # optional: HTML report (needs Quarto)
```
Run from the project root so the `R/` and `figures/` relative paths resolve.

## Where things live
- `R/simulate_data.R` — synthetic inputs (population, rates, comorbidities).
- `R/engine_naive.R` — slow reference engine (string state).
- `R/engine_fast.R` — vectorized engine (logical matrix). The two must stay
  equivalent: after any change to one, re-run `scripts/01_benchmark.R` and
  confirm `max |naive - fast| = 0`.
- `R/psa.R` — Latin Hypercube Sampling + PSA driver + interval summaries.

## Common edits
- **Bigger/smaller run**: change `n` and `n_cycles` in the `make_inputs(...)`
  calls inside the scripts. The naive engine is O(n · cycles · diseases); keep
  `n` modest (≤ 5000) when the naive engine is involved.
- **More/fewer PSA draws**: change `N_DRAWS` in `scripts/02_psa.R`.
- **Add a disease**: add it to `DISEASES` in `R/simulate_data.R` and give it a
  base rate in `build_rates()`. Both engines pick it up automatically.
- **Add a comorbidity**: add a row to `build_comorbidity()` (risk_factor →
  outcome, rr). Re-run the benchmark to confirm equivalence still holds.

## Guardrails
- Keep the transition rule identical in both engines:
  `p = 1 - exp(-(rate * exposure * comorbidity))`.
- Both engines must consume the RNG in the same order to stay bit-identical;
  if you change draw order in one, change it in the other.
- Figures are committed so the README renders on GitHub — regenerate them after
  changing any model code.
