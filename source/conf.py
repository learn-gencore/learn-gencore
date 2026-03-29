project = 'Core BioInformatics  NYU Abu Dhabi'
copyright = '2026, Core Bioinformatics | NYU Abu Dhabi'
author = 'Core Bioinformatics NYU Abu Dhabi'
release = '1.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx_panels']
templates_path = ['_templates']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
import sphinx_rtd_theme
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_show_sourcelink = False
html_show_sphinx = False

# custom.css is inside one of the html_static_path folders (e.g. _static)
html_css_files = ["custom.css"]
html_favicon = "img/favicon.ico"