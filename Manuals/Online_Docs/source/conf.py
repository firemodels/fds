# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'FDS online docs'
copyright = 'https://www.nist.gov/oism/copyrights'
author = 'Jaskaran Gill'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.extlinks',
    'sphinx_rtd_theme'
]

templates_path = ['_templates']
exclude_patterns = []

# Make sure the target is unique
autosectionlabel_prefix_document = True

# Enable numref
numfig = True
math_numfig = True
numfig_secnum_depth = 1

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html #options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = []
html_css_files = [
    'custom.css',
    'css/style.css' 
]

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.

html_logo = 'FDS-logos/FDS-logos_white_cropped.png'
html_favicon = 'FDS-logos/FDS-logos_black_small_square.png'


