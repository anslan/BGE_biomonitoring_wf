# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os

# Tell Jinja2 templates the build is running on Read the Docs
if os.environ.get("READTHEDOCS", "") == "True":
    if "html_context" not in globals():
        html_context = {}
    html_context["READTHEDOCS"] = True

# -- Project information -----------------------------------------------------

project = 'bioscanflow'
copyright = '2026, Sten Anslan'
author = 'Sten Anslan et al.'

# The full version, including alpha/beta/rc tags
version = '1.0b'
release = version

# -- General configuration ---------------------------------------------------
extensions = ['sphinx_design', 'sphinx_copybutton']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build', '_local_docs']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'sphinx_rtd_theme'

html_js_files = ["js/downloadbutton.js"]
html_css_files = ["css/custom.css"]

html_theme_options = {
    'sidebarwidth': 300,
    'collapse_navigation': True
}

sphinx_tabs_disable_tab_closing = True
html_allow_raw_html = True

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

def setup(app):
    app.add_css_file('custom.css')

