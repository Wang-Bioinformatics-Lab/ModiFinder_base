import unittest
import modifinder as mf
from modifinder.convert import (
    to_compound
)

from modifinder.samples import (
    caffeine as caffeine_data
)
import modifinder.tests.utils as utils

class TestConvert(unittest.TestCase):
    def test_to_compound(self):
        compounds = {}
        # compound from accession
        accession = caffeine_data.accession
        compound_from_accession = to_compound(accession)
        self.assertFalse(compound_from_accession is None)
        compounds["from_accession"] = compound_from_accession

        usi = caffeine_data.usi
        compound_from_usi = to_compound(usi)
        self.assertFalse(compound_from_usi is None)
        compounds["from_usi"] = compound_from_usi


        # compound from dict
        data = caffeine_data.data
        compound_from_dict = to_compound(data)
        self.assertFalse(compound_from_dict is None)
        compounds["from_dict"] = compound_from_dict


        # compound from Compound
        compound_from_compound = to_compound(caffeine_data.compound)
        self.assertNotEqual(compound_from_compound, compound_from_accession)
        compounds["from_compound"] = compound_from_compound


        # check if the compound from accession, usi and dict are the same
        for key in ["adduct_mass", "is_known", "name"]:
            self.assertEqual(getattr(compound_from_accession, key), getattr(compound_from_usi, key))
            self.assertEqual(getattr(compound_from_accession, key), getattr(compound_from_dict, key))
            self.assertEqual(getattr(compound_from_accession, key), getattr(compound_from_compound, key))
        
        # check the structures to be the same
        for compound in compounds.values():
            self.assertTrue(compound.structure is not None)
            self.assertTrue(compound.structure.HasSubstructMatch(compound_from_accession.structure))
            self.assertTrue(compound_from_accession.structure.HasSubstructMatch(compound.structure))
        
        # check the spectra to be the same
        for compound in compounds.values():
            self.assertTrue(compound.spectrum is not None)
            self.assertEqual(compound.spectrum.mz, compound_from_accession.spectrum.mz)
            self.assertEqual(compound.spectrum.intensity, compound_from_accession.spectrum.intensity)
            self.assertEqual(compound.spectrum.precursor_mz, compound_from_accession.spectrum.precursor_mz)
            self.assertEqual(compound.spectrum.precursor_charge, compound_from_accession.spectrum.precursor_charge)
            self.assertEqual(compound.spectrum.adduct, compound_from_accession.spectrum.adduct)
    

    def test_to_compound_use_object(self):
        base_compound = caffeine_data.compound
        base_compound.test_arg = "base"
        base_compound_copy = mf.Compound(base_compound)

        self.assertNotEqual(base_compound, base_compound_copy)
        self.assertTrue(utils.compare_compounds(base_compound, base_compound_copy))
        self.assertEqual(base_compound.test_arg, base_compound_copy.test_arg)

        accession = caffeine_data.accession
        compound_from_accession = to_compound(accession, use_object=base_compound)
        compound_from_accession.test_arg = "accession"
        self.assertFalse(compound_from_accession is None)
        self.assertEqual(base_compound, compound_from_accession)
        self.assertEqual(base_compound.test_arg, "accession")
        self.assertEqual(base_compound_copy.test_arg, "base")

        usi = caffeine_data.usi
        compound_from_usi = to_compound(usi, use_object=base_compound)
        compound_from_usi.test_arg = "usi"
        self.assertFalse(compound_from_usi is None)
        self.assertEqual(base_compound, compound_from_usi)
        self.assertEqual(base_compound.test_arg, "usi")
        self.assertEqual(base_compound_copy.test_arg, "base")

        # compound from dict
        data = caffeine_data.data
        compound_from_dict = to_compound(data, use_object=base_compound)
        compound_from_dict.test_arg = "dict"
        self.assertFalse(compound_from_dict is None)
        self.assertEqual(base_compound, compound_from_dict)
        self.assertEqual(base_compound.test_arg, "dict")
        self.assertEqual(base_compound_copy.test_arg, "base")

        # compound from Compound
        compound_from_compound = to_compound(caffeine_data.compound, use_object=base_compound)
        compound_from_compound.test_arg = "compound"
        self.assertEqual(base_compound, compound_from_compound)
        self.assertEqual(base_compound.test_arg, "compound")
        self.assertEqual(base_compound_copy.test_arg, "base")
        self.assertEqual(compound_from_compound, base_compound)
            
            


        # def test_to_compound_from_Compound(self):
        #     compound = mf.Compound()
        #     result = to_compound(compound)
        #     self.assertEqual(result, compound)

        # def test_to_compound_from_USI(self):
        #     usi = "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00005435812"
        #     result = to_compound(usi)
        #     self.assertEqual(result.usi, usi)
        
        # def test_to_compound_from_dict(self):
        #     data = {
        #         "peaks": [{"mz": 100.0, "intensity": 200.0}, {"mz": 150.0, "intensity": 300.0}],
        #         "Precursor_MZ": 500.0,
        #         "Charge": 2
        #     }
        #     result = to_compound(data)
        #     self.assertEqual(result.peaks, data["peaks"])
        #     self.assertEqual(result.precursor_mz, data["Precursor_MZ"])
        #     self.assertEqual(result.charge, data["Charge"])
        
        # def test_to_compound_from_invalid_data(self):
        #     with self.assertRaises(mf.ModiFinderError):
        #         to_compound(123)
        