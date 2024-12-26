import unittest
import json
import modifinder as mf
from modifinder.tests import utils as test_utils
from modifinder import convert as convert
from modifinder.utilities.network import get_matched_peaks
from modifinder.engines.alignment.CosineAlignmentEngine import CosineAlignmentEngine
from modifinder.samples import (
    caffeine as caffeine_data,
    theophylline as theophylline_data
)

class TestCosineAlignment(unittest.TestCase):
    def test_single_alignment(self):
        spectrum1 = convert.to_spectrum(caffeine_data.usi, normalize_peaks=True)
        spectrum2 = convert.to_spectrum(theophylline_data.usi, normalize_peaks=True)
        
        cosine_engine = CosineAlignmentEngine()
        edge_detail = cosine_engine.single_align(spectrum1, spectrum2)
        
        network_match = get_matched_peaks(caffeine_data.usi, theophylline_data.usi)
        
        self.assertAlmostEqual(edge_detail.match_score, network_match["cosine"], places=4)
        self.assertEqual(len(edge_detail.matches), network_match["n_peak_matches"])
    
    def test_align(self):
        modifinder = mf.ModiFinder(caffeine_data.compound, theophylline_data.compound)
        cosine_engine = CosineAlignmentEngine()
        cosine_engine.align(modifinder.network)
        
        if caffeine_data.compound.spectrum.precursor_mz < theophylline_data.compound.spectrum.precursor_mz:
            edge_data =  modifinder.network.get_edge_data(caffeine_data.compound.id, theophylline_data.compound.id)
            network_match = get_matched_peaks(caffeine_data.usi, theophylline_data.usi)
        else:
            edge_data =  modifinder.network.get_edge_data(theophylline_data.compound.id, caffeine_data.compound.id)
            network_match = get_matched_peaks(theophylline_data.usi, caffeine_data.usi)
        
        self.assertAlmostEqual(edge_data["edgedetail"].match_score, network_match["cosine"], places=4)
        self.assertEqual(len(edge_data["edgedetail"].matches), network_match["n_peak_matches"])
        
        
        