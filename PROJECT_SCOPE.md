# Project Scope (Active Production Path)

This repository's **primary purpose** is to produce influenza burden estimates for the
RespiCompass 2024 scenario round using the current ECDC modelling workflow.

## Active production path

The active, supported production path is:

1. `code/00_main.R` (entrypoint)
2. `code/02_settings/settings_version0.R` (scenario and runtime settings)
3. `code/01_main_supporting/load_flu_data.R` (data loading and caching)
4. `code/01_main_supporting/run_flu_models.R` (model orchestration)
5. `code/01_main_supporting/model_SIR_multiseason.R` (SIR multi-season fit/projection)
6. `code/01_main_supporting/process_and_save.R` (post-processing and outputs)

## In scope

- RespiCompass 2024 scenario projections.
- Country-level burden estimation with uncertainty under configured scenarios.
- Operational reproducibility of the above workflow.

## Out of scope (for production)

- General multi-pathogen framework design.
- Experimental or sandbox analyses not required for the RespiCompass 2024 deliverable.
- Legacy modelling scripts that are not invoked by `code/00_main.R`.

## Legacy and exploratory code

Existing legacy and exploratory code may remain in the repository for reference,
but is not considered production unless explicitly linked into the active path above.
