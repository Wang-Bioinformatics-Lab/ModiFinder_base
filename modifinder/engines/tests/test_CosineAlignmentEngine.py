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
        spectrum1 = convert.to_spectrum(caffeine_data.spectrum)
        spectrum2 = convert.to_spectrum(theophylline_data.spectrum)

        cosine_engine = CosineAlignmentEngine()
        edge_detail = cosine_engine.align_single(
            spectrum1,
            spectrum2,
            caffeine_data.compound.id,
            theophylline_data.compound.id
            )
        
        network_match = get_matched_peaks(
            caffeine_data.compound.id,
            theophylline_data.compound.id,
        )
        
        self.assertAlmostEqual(edge_detail.match_score, network_match["cosine"], places=4)
        self.assertEqual(len(edge_detail.matches), network_match["n_peak_matches"])
    
    @unittest.expectedFailure
    def test_align(self):
        """This test is an expected failure, precursor m/z filtering happens automatically when a ModiFinder Compound
        is created, altering the cosine similarity.        
        """
        modifinder = mf.ModiFinder(
            caffeine_data.compound,
            theophylline_data.compound,
            ppm_tolerance = 1000    # Since resolver only does m/z, we should set ppm to very high
        )
        cosine_engine = CosineAlignmentEngine()

        cosine_engine.align(
            modifinder.network,
            mz_tolerance=0.1,
            align_all=False
        )
        
        if caffeine_data.compound.spectrum.precursor_mz < theophylline_data.compound.spectrum.precursor_mz:
            edge_data =  modifinder.network.get_edge_data(caffeine_data.compound.id, theophylline_data.compound.id)
            network_match = get_matched_peaks(caffeine_data.compound.id, theophylline_data.compound.id)
        else:
            edge_data =  modifinder.network.get_edge_data(theophylline_data.compound.id, caffeine_data.compound.id)
            network_match = get_matched_peaks(theophylline_data.compound.id, caffeine_data.compound.id)
        
        self.assertAlmostEqual(edge_data["edgeDetail"].match_score, network_match["cosine"], places=4)
        self.assertEqual(len(edge_data["edgeDetail"].matches), network_match["n_peak_matches"])
        
        
    def test_align_manual(self):
        """This test is an expected failure, precursor m/z filtering happens automatically when a ModiFinder Compound
        is created, altering the cosine similarity.        
        """
        modifinder = mf.ModiFinder(
            caffeine_data.compound,
            theophylline_data.compound,
            ppm_tolerance = 1000    # Since resolver only does m/z, we should set ppm to very high
        )
        cosine_engine = CosineAlignmentEngine()

        cosine_engine.align(
            modifinder.network,
            mz_tolerance=0.1,
            align_all=False
        )
        
        if caffeine_data.compound.spectrum.precursor_mz < theophylline_data.compound.spectrum.precursor_mz:
            edge_data =  modifinder.network.get_edge_data(caffeine_data.compound.id, theophylline_data.compound.id)
        else:
            edge_data =  modifinder.network.get_edge_data(theophylline_data.compound.id, caffeine_data.compound.id)
        
        self.assertAlmostEqual(edge_data["edgeDetail"].match_score, 0.9868885421288105, places=4)
        self.assertEqual(len(edge_data["edgeDetail"].matches), 1)