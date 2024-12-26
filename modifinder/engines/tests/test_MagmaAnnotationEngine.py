import unittest
import json
import modifinder as mf
from modifinder.tests import utils as test_utils
from modifinder import convert as convert
from modifinder.utilities.network import get_matched_peaks
from modifinder.engines.annotation.MAGMaAnnotationEngine import MAGMaAnnotationEngine
from modifinder.samples import (
    caffeine as caffeine_data,
    theophylline as theophylline_data
)
from modifinder.engines.annotation.magma import (
    fragmentation_py,
    rdkit_engine as rdkit_engine,
)

class TestMagmaAnnotation(unittest.TestCase):
    def test_single_annotation(self):
        compound = caffeine_data.compound
        breaks = 4
        magma_engine = MAGMaAnnotationEngine(breaks=breaks)
        magma_engine.annotate_single(compound, modify_compound=True, ppm=45)
        self.assertTrue(compound.peak_fragments_map is not None)    
            
                
                
        
    