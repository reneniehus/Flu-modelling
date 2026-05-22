# Production Contract (RespiCompass 2024 workflow)

## What "contract" means here

In this context, a **contract** is an explicit set of expectations that must stay
true while refactoring:

- which files define the production workflow,
- what inputs they require,
- what outputs they produce,
- and what behaviors must not change without deliberate review.

This makes architecture changes safer by giving us a stable baseline to test against.

## Contract: execution flow

The production flow is:

1. Source setup and support scripts from `code/00_main.R`.
2. Build settings from `code/02_settings/settings_version0.R`.
3. Load data with `load_flu_data(params, regenerate = F, new_from_online = F)`.
4. Run model orchestration with `run_flu_models(params, data)`.
5. Postprocess and outputs with `process_and_save(...)`.

## Contract: model selection and settings

- Primary model: `SIR_simple_multi_season`.
- Current settings include practical fast-path controls such as:
  - `params$load_earlyfit` (load cached fit when available),
  - `params$rapid_stan_fit` (lighter fitting behavior when fitting is performed),
  - `params$run_countries` (country scope restriction).

## Contract: data behavior

- Local/cached data loading behavior should remain available and predictable.
- Data loading should respect caller-specified `new_from_online` and `regenerate`
  behavior in the active pipeline.

## Contract: outputs

The post-processing stage must continue to provide:

- submission-ready scenario burden outputs,
- parameter summary outputs,
- core summary objects required for reporting.

## Change policy

Any change that alters this contract should be:

1. intentional,
2. documented,
3. reviewed with project-purpose rationale,
4. validated with a lightweight check before broader refactoring.
