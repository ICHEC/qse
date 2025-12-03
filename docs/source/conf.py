# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

import tomli

sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

with open("../../pyproject.toml", "rb") as f:
    toml_dict = tomli.load(f)

project = toml_dict["project"]["name"]
copyright = "..."
author = "ICHEC"
release = toml_dict["project"]["version"]

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "nbsphinx",
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_book_theme",
    "sphinx_jupyterbook_latex",
    "sphinx_togglebutton",
    "sphinx_copybutton",
]
nbsphinx_execute = 'never'
templates_path = ["_templates"]
exclude_patterns = []


# Napoleon settings
napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_use_param = False
napoleon_use_ivar = True
napoleon_use_rtype = False

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_book_theme"
html_static_path = ["_static"]
html_theme_options = {
    "repository_url": "https://github.com/ICHEC/qse",
    "repository_branch": "main",
    "use_repository_button": True,
    "use_edit_page_button": False,
    "use_issues_button": True,
}
html_logo = "logo.svg"
