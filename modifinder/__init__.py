"""
ModiFinder
==========

ModiFinder is a Python package for the identification of modifications in mass spectrometry data.
"""

__version__ = "1.4.62beta"

import modifinder.convert as convert
from modifinder.classes import Compound, Spectrum, ModiFinder, EdgeDetail, MatchType, StructureMeta
from modifinder.convert import (
    to_compound,
    to_spectrum,
    compound_to_dict,
    spectrum_to_dict,
)
from modifinder.exceptions import (
    ModiFinderError,
    ModiFinderException,
    ModiFinderNetworkError,
    ModiFinderNotImplementedError,
    ModiFinderNotSolvableError,
)
from modifinder.engines.alignment.CosineAlignmentEngine import CosineAlignmentEngine
from modifinder.engines.annotation.MAGMaAnnotationEngine import MAGMaAnnotationEngine
from modifinder.engines.evaluation.BasicEvaluationEngine import BasicEvaluationEngine

__all__ = [
    "__version__",
    "Compound",
    "Spectrum",
    "ModiFinder",
    "EdgeDetail",
    "MatchType",
    "StructureMeta",
    "convert",
    "to_compound",
    "to_spectrum",
    "compound_to_dict",
    "spectrum_to_dict",
    "ModiFinderException",
    "ModiFinderError",
    "ModiFinderNetworkError",
    "ModiFinderNotImplementedError",
    "ModiFinderNotSolvableError",
    "CosineAlignmentEngine",
    "MAGMaAnnotationEngine",
    "BasicEvaluationEngine",
]
