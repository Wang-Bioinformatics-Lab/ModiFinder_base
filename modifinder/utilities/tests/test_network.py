import unittest
import json
import modifinder as mf
from modifinder.utilities.network import (
    get_data,
    get_matched_peaks
)

from modifinder.samples import caffeine as caffeine_data

from modifinder.exceptions import ModiFinderNetworkError

class TestConvert(unittest.TestCase):
    def test_get_data(self):
        # test with a valid USI
        usi = caffeine_data.usi
        data_usi = get_data(usi)
        self.assertTrue(data_usi)

        # test with a valid Accession
        accession = caffeine_data.accession
        data_accession = get_data(accession)
        self.assertTrue(data_accession)

        # test both data are the same
        self.assertEqual(data_usi, data_accession)

        # test with an invalid USI
        usi = "mzspec:GNPS:UNKNOWN:accession:UNKNOWN"
        self.assertRaises(ModiFinderNetworkError, get_data, usi)
    
    def test_partial_data(self):
        usi = "mzspec:GNPS2:TASK-b6e775882a194d0185ff8ddf4fd048a6-inputspectra/TUEBINGEN-NATURAL-PRODUCT-COLLECTION-unknowns.mgf:scan:2471"
        data = get_data(usi)
        self.assertTrue(data)
        self.assertTrue("peaks" in data)
        self.assertTrue("precursor_charge" in data)
        self.assertTrue("precursor_mz" in data)



