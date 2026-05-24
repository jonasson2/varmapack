# Configuration file for the Sphinx documentation builder.

project = "Varmapack"
copyright = "2026, Kristján Jónasson"
author = "Kristján Jónasson"

import pathlib
import tomllib

pyproject = pathlib.Path(__file__).parent.parent / "pyproject.toml"
data = tomllib.loads(pyproject.read_text(encoding="utf-8"))
release = data["project"]["version"]

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "numpydoc",
]

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
