# Configuration file for the Sphinx documentation builder.

project = "Varmapack"
copyright = "2026, Kristján Jónasson"
author = "Kristján Jónasson"

import os
import pathlib
import sys
import tomllib

docs_dir = pathlib.Path(__file__).parent
python_dir = pathlib.Path(__file__).parent.parent
if os.environ.get("READTHEDOCS") == "True":
    sys.path.insert(0, str(docs_dir / "_stubs"))
else:
    sys.path.insert(0, str(python_dir))
pyproject = python_dir / "pyproject.toml"
data = tomllib.loads(pyproject.read_text(encoding="utf-8"))
release = data["project"]["version"]

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "numpydoc",
]

mathjax3_config = {
    "tex": {
        "macros": {
            "eps": r"\varepsilon",
            "Cov": r"\operatorname{Cov}",
            "Var": r"\operatorname{Var}",
        }
    }
}

autosummary_generate = True
autosummary_imported_members = True
numpydoc_show_class_members = False
add_module_names = False
autodoc_typehints = "none"
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "furo"
html_theme_options = {
    "sidebar_hide_name": True,
}
html_static_path = ["_static"]
html_css_files = ["custom.css"]
pygments_style = "default"
pygments_dark_style = "monokai"
