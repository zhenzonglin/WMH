# WMH analysis project

- Work in the existing clone. All statistical code is Python; use `uv sync --frozen` and the locked Python 3.11 environment.
- Patient data are external and read-only. Never commit raw SAS, patient CSV, images, ID maps, local configurations, or outputs. Public SAS fixtures under `tests/fixtures/sas` are the only bundled SAS data.
- The workstation configuration is `config/workstation.local.yml`. Start with `wmh-hcy audit`; accept the actual SuStaIn path with `configure --sustain-dir`; then run through `prepare` before fitting models.
- Follow H1 through H4 and `docs/statistical_analysis_plan.md`. Use only the source fields in `fields.py`; do not invent unavailable variables, infer a center, change the endpoint, silently drop required covariates, or substitute synthetic data for real input.
- Review `outputs/real/audit/` and `prepared/flow.csv` before a real statistical run. Preserve exact IDs, QC decisions, and reasons for exclusions.
- Keep method definitions and CLI documentation consistent. After changes, run `uv run ruff check src tests` and `uv run pytest -q`; add targeted verification for changed statistical or data-handling behavior.
- Report software validation separately from patient-level results. Synthetic effects and counts are not research findings.
