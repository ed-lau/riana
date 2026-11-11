"""Documentation that is generated from the code must be current."""

import runpy
from pathlib import Path

GENERATOR = Path(__file__).resolve().parents[1] / "riana_website" / "_gen_reference.py"


def test_cli_reference_page_is_current():
    """riana_website/reference.qmd is generated from the CLI's help text; regenerate it
    with `python riana_website/_gen_reference.py` after changing an option."""
    gen = runpy.run_path(str(GENERATOR))
    committed = Path(gen["PAGE"]).read_text()
    assert committed == gen["build_reference"](), (
        "riana_website/reference.qmd is out of date with the CLI — run "
        "`python riana_website/_gen_reference.py` and commit the result.")
