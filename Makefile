.PHONY: install test quick production track-plan

install:
	python -m pip install -e .

test:
	python -m unittest discover -s tests -v

quick:
	python -m rfq_pipeline --config configs/isac2_quick.toml all

production:
	python -m rfq_pipeline --config configs/isac2_production.toml all

track-plan:
	python -m rfq_pipeline --config configs/isac2_production.toml track --dry-run

