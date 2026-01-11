Tutorial
========

ModiFinder is a powerful toolkit for mass spectrometry data analysis. While it excels at finding structural modification sites, its utilities are valuable for any MS data analysis workflow. This tutorial covers both the core ModiFinder functionality and the visualization and data processing tools available in the package.


The following tutorials guide you through different aspects of ModiFinder:

**Basics** - Start here to learn ModiFinder's core functionality: comparing MS/MS spectra to predict modification sites in molecules. You'll learn how to create Spectrum and Compound objects, use USI identifiers, and run site localization analysis.

**Customization** - Dive deeper into customizing ModiFinder's behavior. Learn how to adjust scoring parameters, tolerance settings, and prediction strategies to optimize results for your specific compounds and experimental conditions.

**Utilities** - Discover ModiFinder's essential data processing tools: fetching data from GNPS, creating consensus spectra from replicates, reading/writing MGF files, calculating molecular properties, and measuring structural similarity between molecules.

**Visualization** - Learn how to visualize your results effectively: draw molecules with highlighted modification sites, create annotated spectra plots, generate mirror plots for spectral comparison, and produce publication-ready figures.

.. toctree::
    :maxdepth: 2
    
    tutorials/basics.rst
    tutorials/customization.rst
    tutorials/utilities.rst
    tutorials/visualization.rst
