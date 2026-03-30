"""Launch the Streamlit dashboard: python -m dashboard"""

import sys
from pathlib import Path

from streamlit.web.cli import main as st_main

sys.argv = ["streamlit", "run", str(Path(__file__).with_name("app.py"))]
st_main()
