import os
import sys

import tomli

sys.path.insert(0, os.path.abspath("../"))
sys.path.insert(0, os.path.abspath("../../"))

with open("../../pyproject.toml", "rb") as f:
    toml_dict = tomli.load(f)

project = toml_dict["project"]["name"]
copyright = "..."
author = "ICHEC"
release = toml_dict["project"]["version"]

# autosummary_generate = True
bibtex_bibfiles = []
comments_config = {"hypothesis": False, "utterances": False}
exclude_patterns = ["**.ipynb_checkpoints", ".DS_Store", "Thumbs.db", "_build"]

extensions = [
    "sphinx_togglebutton",
    "sphinx_copybutton",
    "myst_nb",
    "jupyter_book",
    "sphinx_thebe",
    "sphinx_comments",
    "sphinx_external_toc",
    "sphinx.ext.intersphinx",
    "sphinx_design",
    "sphinx_book_theme",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinxcontrib.bibtex",
    "sphinx_jupyterbook_latex",
    "sphinx_multitoc_numbering",
    "autoapi.extension",
]
autoapi_dirs = ["../../qse"]  # for sphinx-autoapi
autoapi_options = [  # https://sphinx-autoapi.readthedocs.io/en/latest/reference/config.html#confval-autoapi_options
    "members",
    "show-module-summary",
    "show-inheritance",
]

external_toc_exclude_missing = False
external_toc_path = "_toc.yml"
html_baseurl = ""
html_favicon = ""
html_sourcelink_suffix = ""
html_theme = "sphinx_book_theme"
html_theme_options = {
    "search_bar_text": "Search this book...",
    "launch_buttons": {
        "notebook_interface": "classic",
        "binderhub_url": "",
        "jupyterhub_url": "",
        "thebe": False,
        "colab_url": "",
        "deepnote_url": "",
    },
    "path_to_docs": "docs",
    "repository_url": "https://github.com/ICHEC/qse",
    "repository_branch": "main",
    "extra_footer": "",
    "home_page_in_toc": True,
    "announcement": "",
    "analytics": {
        "google_analytics_id": "",
        "plausible_analytics_domain": "",
        "plausible_analytics_url": "https://plausible.io/js/script.js",
    },
    "use_repository_button": True,
    "use_edit_page_button": False,
    "use_issues_button": True,
}
html_title = "QSE"
latex_engine = "pdflatex"
myst_enable_extensions = [
    "colon_fence",
    "dollarmath",
    "linkify",
    "substitution",
    "tasklist",
]
myst_url_schemes = ["mailto", "http", "https"]
nb_execution_allow_errors = False
nb_execution_cache_path = ""
nb_execution_excludepatterns = []
nb_execution_in_temp = False
nb_execution_mode = "force"
nb_execution_timeout = 30
nb_output_stderr = "show"
numfig = True
pygments_style = "sphinx"
suppress_warnings = ["myst.domains"]
use_jupyterbook_latex = True
use_multitoc_numbering = True

# napoleon settings
napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_use_param = False
napoleon_use_ivar = True
