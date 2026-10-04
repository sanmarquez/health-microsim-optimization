# Project notes for Claude Code

This is a **synthetic, self-contained showcase** of a health micro-simulation:
a slow reference engine vs a vectorized engine (same results, ~100–300× faster),
plus an LHS/PSA uncertainty layer. Base R only for the scripts; Quarto only for
the report. No real data is used anywhere.

## How to run
```bash
Rscript scripts/01_benchmark.R   # speed + equivalence
Rscript scripts/02_psa.R         # LHS/PSA intervals
quarto render report.qmd         # optional HTML report
```
Always run from the project root.

## The one invariant that matters
`engine_naive.R` and `engine_fast.R` model the **same thing** and must produce
**identical** output (`max |naive - fast| = 0`). The fast engine is the point of
the repo; the naive one is the reference that proves it is correct. After editing
either engine, re-run `scripts/01_benchmark.R` and check the equivalence line.

## Style
- Keep it dependency-free (base R) so anyone can run it.
- Keep the two engines readable; this repo is meant to be *read*, not just run.
- Regenerate `figures/` after any model change (the README embeds them).

A `run-pipeline` skill in `.claude/skills/` has the detailed steps.
