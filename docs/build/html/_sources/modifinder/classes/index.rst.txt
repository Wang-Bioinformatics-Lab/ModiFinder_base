Classes
=======

The following classes are defined in the modifinder package:

* Compound: Represents a compound, contains the spectrum, structure, annotation of the peaks, ...
* Spectrum: Represents a spectrum, contains the mz values, intensities, adduct, precursor_mz, ...
* EdgeDetail: Represents an alignment, contains the alignment score, matched peaks and the type of the match
* ModiFinder: The main class of the package, using the known analogs, tries to provide structural information about the unknown compound.

.. toctree::
   :maxdepth: 1

   classes.Spectrum
   classes.Compound
   classes.EdgeDetail
   classes.ModiFinder