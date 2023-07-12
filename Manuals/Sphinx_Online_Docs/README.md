# Welcome to FDS Online Documentation

To compile the online docs, please follow the steps below.

The basic process is to convert the LaTeX manuals to restructured text (.rst) files using [Pandoc](https://pandoc.org/), do a bit of cleanup, and then generate the html using [Sphinx](https://www.sphinx-doc.org/).

Prerequisites:

* Install Pandoc
* Install Sphinx v5.0.2 (Anaconda latest)
* Install [Read the Docs](https://sphinx-rtd-theme.readthedocs.io/en/stable/) theme
* Populate `SCRIPT_FIGURES` directories in `fds/Manuals/*/`

Steps to compile:

1. Run Pandoc script (under construction)
2. Run Sphinx `make html`
