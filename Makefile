.PHONY: run dashboard fmt lint check install install-dev install-notebooks

# Install core dependencies
install:
	pip install -e .

# Install with notebook dependencies (matplotlib, pandas, dlib)
install-notebooks:
	pip install -e ".[notebooks]"

# Install with dev tools (ruff)
install-dev:
	pip install -e ".[dev]"

# Run the CLI simulation
run:
	python -m app

# Launch the interactive Streamlit dashboard
dashboard:
	streamlit run dashboard/app.py

# Format code with ruff
fmt:
	ruff format app/ dashboard/

# Lint code with ruff
lint:
	ruff check app/ dashboard/

# Lint and auto-fix
fix:
	ruff check --fix app/ dashboard/

# Format + lint (pre-commit check)
check: fmt lint
