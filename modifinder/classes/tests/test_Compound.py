import unittest
import json
import modifinder as mf
from modifinder.tests import utils as test_utils
from modifinder import Compound
from modifinder.samples import (
    caffeine as caffeine_data
)

class TestConvert(unittest.TestCase):
    def test_create_compound(self):
        peaks = json.loads(caffeine_data.data["peaks_json"])
        precursor_mz = caffeine_data.data["Precursor_MZ"]
        precursor_charge = caffeine_data.data["Charge"]
        adduct = caffeine_data.data["Adduct"]
        smiles = caffeine_data.data["Smiles"]
        name = caffeine_data.data["Compound_Name"]
        compound = Compound(id="CCMSLIB00005435812", peaks=peaks, precursor_mz=precursor_mz, precursor_charge=precursor_charge, adduct=adduct, smiles=smiles, name=name)
        # print("done making the compound with all the data")
        # print("=====================================================")
        compound_from_accession = Compound(caffeine_data.accession)
        # print("done making the compound from accession")
        # print("=====================================================")

        compound_from_usi = Compound(caffeine_data.usi)

        data_dict = {
            "id": "CCMSLIB00005435812",
            "peaks": peaks,
            "precursor_mz": precursor_mz,
            "charge": precursor_charge,
            "adduct": adduct,
            "smiles": smiles,
            "name": name
        }
        compound_from_dict = Compound(data_dict)

        self.assertTrue(test_utils.compare_compounds(compound, compound_from_accession))
        self.assertTrue(test_utils.compare_compounds(compound, compound_from_usi))
        self.assertTrue(test_utils.compare_compounds(compound, compound_from_dict))
    
    def test_clear(self):
        peaks = json.loads(caffeine_data.data["peaks_json"])
        precursor_mz = caffeine_data.data["Precursor_MZ"]
        precursor_charge = caffeine_data.data["Charge"]
        adduct = caffeine_data.data["Adduct"]
        smiles = caffeine_data.data["Smiles"]
        name = caffeine_data.data["Compound_Name"]
        data_dict = {
            "id": "CCMSLIB00005435812",
            "peaks": peaks,
            "precursor_mz": precursor_mz,
            "charge": precursor_charge,
            "adduct": adduct,
            "smiles": smiles,
            "name": name
        }
        
        compound_from_accession = Compound(caffeine_data.accession)
        # check for additional argument create_time to exist here
        self.assertTrue(hasattr(compound_from_accession, "create_time"))
        compound_from_accession.clear()
        compound_from_accession.update(**data_dict)
        # check for additional argument create_time to not exist here
        self.assertFalse(hasattr(compound_from_accession, "create_time"))

