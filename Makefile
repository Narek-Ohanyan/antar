.PHONY: install test figures audit
install:
	pip install -e ".[dev]"
test:
	python -m pytest
figures:
	python scripts/make_figures.py --out docs/figures
audit:
	python scripts/audit_v1.py --parquet ../armenia-reforestation-predictive-framework-main/data/Armenia_ML_Training_Data.parquet
