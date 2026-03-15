"""Sphinx configuration for gravica documentation."""

project = "gravica"
author = "gravica contributors"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx_autodoc_typehints",
]

html_theme = "furo"

autodoc_member_order = "bysource"
