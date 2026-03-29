.PHONY: run dashboard fmt lint check install install-dev install-all

# Install core dependencies
install:
	pip install -e .

# Install with dev tools (ruff)
install-dev:
	pip install -e ".[dev]"

# Install everything (core + dashboard + dev)
install-all:
	pip install -e ".[dashboard,dev]"

# Run the CLI simulation (batch mode, saves PNGs to output/)
run:
	python -m app

# Launch the interactive Streamlit dashboard
dashboard:
	streamlit run app/dashboard.py

# Format code with ruff
fmt:
	ruff format app/

# Lint code with ruff
lint:
	ruff check app/

# Lint and auto-fix
fix:
	ruff check --fix app/

# Format + lint (pre-commit check)
check: fmt lint
