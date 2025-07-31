# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import pathlib
import sys
import logging

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

logger.info(f"python exec: {sys.executable}")
logger.info(f"sys.path: {sys.path}")
root = pathlib.Path(__file__).parent.parent.absolute()
os.environ["PYTHONPATH"] = str(root)
sys.path.insert(0, str(root))

project = "ROMS"
copyright = "2025, UCLA, [C]Worthy"
author = "ROMS developers"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinxcontrib.bibtex",
]

templates_path = ["_templates"]
exclude_patterns = []

# Markdown options
myst_enable_extensions = ["dollarmath", "deflist", "colon_fence"] 
myst_heading_anchors = 3  # add anchor links to headers

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_book_theme"
# html_theme = 'alabaster'
html_static_path = ["_static"]

bibtex_bibfiles = ["references.bib"]
bibtex_reference_style = "author_year"

html_theme_options = {
    "repository_url": "https://github.com/CWorthy-ocean/ucla-roms.git",
    "use_repository_button": True,
}
