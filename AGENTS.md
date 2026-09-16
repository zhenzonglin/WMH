# WMH analysis project

- Work in the existing clone. All statistical code is Python 3.11. Support either `conda env create -f environment.yml` or `uv sync --frozen`. Conda uses the version pins exported from `uv.lock` into `requirements-conda.txt`; do not independently upgrade analysis dependencies.
- Patient data are external and read-only. Never commit raw SAS, patient CSV, images, ID maps, local configurations, or outputs. Public SAS fixtures under `tests/fixtures/sas` are the only bundled SAS data.
- The workstation configuration is `config/workstation.local.yml`. Start with `wmh-hcy audit`; accept the actual SuStaIn path with `configure --sustain-dir`; then run through `prepare` before fitting models.
- Follow H1 through H4 and `docs/statistical_analysis_plan.md`. Use only the source fields in `fields.py`; do not invent unavailable variables, infer a center, change the endpoint, silently drop required covariates, or substitute synthetic data for real input.
- Review `outputs/real/audit/` and `prepared/flow.csv` for data availability. Per the researcher's 2026-09-16 revision, manual image review labels do not determine eligibility in any analysis; do not reintroduce a QC gate or relabel records as passed. Preserve exact IDs and reasons for data exclusions.
- Keep method definitions and CLI documentation consistent. In the activated Conda environment run `ruff check src tests` and `python -m pytest -q`; for uv prefix those commands with `uv run`. Add targeted verification for changed statistical or data-handling behavior.
- Report software validation separately from patient-level results. Synthetic effects and counts are not research findings.
