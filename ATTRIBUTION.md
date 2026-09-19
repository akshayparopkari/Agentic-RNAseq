# Attribution

This project builds on the published 3′ Tag-Seq analysis pipeline:

> Paropkari, A. D., Bapat, P. S., Sindi, S. S., & Nobile, C. J. (2023).
> A computational workflow for analysis of 3′ Tag-Seq data.
> *Current Protocols*, 3, e664. https://doi.org/10.1002/cpz1.664

Source repository: https://github.com/akshayparopkari/RNAseq

## Files adapted for this project

- `pipeline/deseq.R` — adapted from the original repo's `deseq.R`. The
  hardcoded reference level (`levels = c("WT", "Mutant")`) is replaced
  with a configurable parameter, chosen by the agent's planning step
  based on the actual `Condition` values present in the sample metadata,
  rather than assumed in advance.
- `pipeline/format_counts_table.py` — vendored unmodified.

## What's untouched

`pipeline.sh`, the original `deseq.R`, `overlap_upsetR.R`, and everything
else in the original repository remains exactly as published there, none
of it is copied into this project beyond the two files listed above. That
repository is co-authored and peer-reviewed work with its own citation,
this project is a separate, personal exploration built on top of it, not
a modification to it.
