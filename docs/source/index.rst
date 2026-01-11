.. ModiFinder documentation master file, created by
   sphinx-quickstart on Fri Oct 25 15:47:52 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ModiFinder's documentation!
======================================

**ModiFinder** is a powerful Python toolkit for mass spectrometry data analysis, specializing in site localization of structural modifications using MS/MS data.

Beyond Modification Finding
----------------------------

While ModiFinder excels at pinpointing where structural modifications occur in molecules, it's also a comprehensive toolkit for:

* **Working with Mass Spectra**: Parse USIs, fetch data from GNPS, and manipulate spectrum objects
* **Molecular Visualization**: Draw molecules, spectra, structural comparisons, and modification heatmaps
* **Data Processing**: Read/write MGF files, normalize spectra, and batch process MS data
* **Compound Management**: Create and manipulate compound objects from SMILES, InChI, or GNPS identifiers
* **Network Analysis**: Retrieve and compare spectra from public repositories

Whether you're analyzing a single modification site or building a large-scale MS data processing pipeline, ModiFinder's utilities can help streamline your workflow.

Quick Start
-----------

Install ModiFinder:

.. code-block:: bash

   pip install modifinder

Get started with a simple example:

.. code-block:: python

   from modifinder import Compound
   from modifinder.utilities import visualizer as viz
   
   # Fetch compound from GNPS
   compound = Compound("CCMSLIB00010113829")
   
   # Draw the structure
   img = viz.draw_molecule(compound.structure)
   
   # Access spectrum data
   print(f"{len(compound.spectrum.mz)} peaks")

Key Features
------------

🎯 **Modification Site Prediction**
   Identify the most likely modification sites in unknown compounds using MS/MS data

🔬 **Spectrum Analysis**
   Parse, normalize, and compare mass spectra from various sources

🎨 **Rich Visualization**
   Create publication-quality figures of molecules, spectra, and structural comparisons

📊 **Data Integration**
   Seamlessly work with GNPS data, MGF files, and USIs

🧪 **Flexible Objects**
   Compound and Spectrum classes that work with your existing data formats

.. note::
   This project is under active development.


Citing
------

.. include:: ../../CITATION.md


License
-------
   
.. include:: ../../LICENSE.md


.. toctree::
   :maxdepth: 4
   :hidden:

   getting_started
   tutorial
   modifinder/index
   

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

