# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# -- Project information -----------------------------------------------------

project = 'bioscanflow'
copyright = '2024, Sten Anslan'
author = 'Sten Anslan et al.'

# The full version, including alpha/beta/rc tags
version = '0.2b'
release = version

# -- General configuration ---------------------------------------------------
# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = []

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', '_local_docs']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'sphinx_rtd_theme'

html_theme_options = {
    'sidebarwidth': 300,
    'collapse_navigation': True
}

sphinx_tabs_disable_tab_closing = True

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# custom css
def setup(app):
    app.add_css_file('text-backgrounds.css')
# usage:
#.. raw:: html

#  <span class="highlight-blue">This text has blue background</span>
