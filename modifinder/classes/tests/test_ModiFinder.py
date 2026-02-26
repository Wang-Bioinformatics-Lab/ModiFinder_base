import unittest
import json
from modifinder.tests import utils as test_utils
from modifinder import ModiFinder
from modifinder.samples import (
    caffeine as caffeine_data,
    theophylline as theophylline_data
)
from modifinder.engines.alignment.CosineAlignmentEngine import CosineAlignmentEngine
from modifinder.engines.annotation.MAGMaAnnotationEngine import MAGMaAnnotationEngine

class TestConvert(unittest.TestCase):
    def test_create_modifinder_usecase1(self):

        def check_network(network, bigger_id, smaller_id):
            self.assertIsNotNone(network)
            self.assertEqual(len(network.nodes), 2)
            self.assertEqual(len(network.edges), 1)
            self.assertTrue(smaller_id in network.nodes)
            self.assertTrue(bigger_id in network.nodes)

            edge = network.edges[(smaller_id, bigger_id)]
            self.assertIsNotNone(edge)
            # self.assertTrue("edgeDetails" in edge)
            # self.assertIsNotNone(edge["edgeDetails"])
            # self.assertEqual(edge["edgeDetails"].number_of_modifications, 1)

        #generate with compound objects
        modifinder = ModiFinder(knownCompound=caffeine_data.compound, unknownCompound=theophylline_data.compound)
        import sys
        print("Nodes in the modifinder.network:", list(modifinder.network.nodes), file=sys.stderr, flush=True)
        print("Edges in the modifinder.network:", list(modifinder.network.edges), file=sys.stderr, flush=True)
        check_network(modifinder.network, caffeine_data.accession, theophylline_data.accession)
    
    def test_generate_probabilities(self):
        modifinder = ModiFinder(knownCompound=caffeine_data.compound, unknownCompound=theophylline_data.compound, ppm = 50)
        print(modifinder.network.edges[(theophylline_data.accession, caffeine_data.accession)]["edgeDetail"])
        print(modifinder.network.nodes[caffeine_data.accession]["compound"].spectrum.peak_fragment_dict)
        print(modifinder.generate_probabilities())


