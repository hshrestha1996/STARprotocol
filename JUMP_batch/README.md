# JUMP batch (`jump -batch-id` / `jump -batch-q`)

Cross-batch peptide identification and TMT quantification for multi-batch
experiments, part of the JUMP Software Suite.

## Layout

```
bin/batch_id.pl              entry point for `jump -batch-id`
bin/batch_q.pl               entry point for `jump -batch-q`
bin/batch/BatchLib.pl        shared subroutines (params, job submission)
bin/batch/jump_f_batch.pl    cross-batch identification rescue
bin/batch/jump_batch_v1.6.pl cross-batch quantification merge
bin/batch/jump_i.pl          QC heatmaps / clustering
bin/batch/R/batch_normalization.R   internal-standard and limma normalisation
lib/                         Perl modules
examples/                    example parameter files and a run log
docs/                        user manual
reference/                   the `jump` driver and `jump_params`, for context
```

`bin/batch_id.pl` resolves `BatchLib.pl` via `$Bin/batch/`, and
`jump_batch_v1.6.pl` resolves the R script via `$Bin/R/`. These relative
positions must be preserved.

## Requirements

- Perl 5.10.1 — `Clone`, `Storable`, `Statistics::R`, `List::Util`,
  `FindBin`, `File::Basename`, `Cwd`, `autodie`
- R — `limma` (Bioconductor)
- A working JUMP installation providing `jump -params`, `-f`, `-q`, `-l`

## Status

See `MISSING.md`. This tree is not yet runnable as-is.
