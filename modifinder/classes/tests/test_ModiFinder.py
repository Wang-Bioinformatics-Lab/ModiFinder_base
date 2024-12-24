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
            # self.assertTrue("edgedetails" in edge)
            # self.assertIsNotNone(edge["edgedetails"])
            # self.assertEqual(edge["edgedetails"].number_of_modifications, 1)

        #generate with compound objects
        modifinder = ModiFinder(knownCompond=caffeine_data.compound, unknownCompound=theophylline_data.compound)
        check_network(modifinder.network, caffeine_data.accession, theophylline_data.accession)


        #generate with usi strings
        modifinder = ModiFinder(knownCompond=caffeine_data.usi, unknownCompound=theophylline_data.usi)
        check_network(modifinder.network, caffeine_data.accession, theophylline_data.accession)

        #generate with accessions
        modifinder = ModiFinder(knownCompond=caffeine_data.accession, unknownCompound=theophylline_data.accession)
        check_network(modifinder.network, caffeine_data.accession, theophylline_data.accession)

        #generate with dictionaries
        modifinder = ModiFinder(knownCompond=caffeine_data.data, unknownCompound=theophylline_data.data)
        check_network(modifinder.network, caffeine_data.accession, theophylline_data.accession)
    
    
    def test_generate_probabilities(self):
        modifinder = ModiFinder(knownCompond=caffeine_data.compound, unknownCompound=theophylline_data.compound, ppm = 50)
        print(modifinder.network.edges[(theophylline_data.accession, caffeine_data.accession)]["edgedetail"])
        print(modifinder.network.nodes[caffeine_data.accession]["compound"].peak_fragments_map)
        print(modifinder.generate_probabilities())


