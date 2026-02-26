"""Test cases using USI-based spectra from GNPS2.

Case 1: Known compound (scan:12430) with SMILES + Modified compound (scan:5538).
Case 2: Known compound (scan:1850) with SMILES + Modified compound (scan:1728).
         This is the typical dashboard use case.

Spectrum data is stored as a JSON fixture in test_files/usi_test_spectra.json
to avoid API calls during testing.
"""

import unittest
import json
import os
import numpy as np

from modifinder import Compound, ModiFinder

FIXTURE_PATH = os.path.join(os.path.dirname(__file__), "test_files", "usi_test_spectra.json")

with open(FIXTURE_PATH) as f:
    _FIXTURES = json.load(f)


def _load_fixture(key: str) -> dict:
    """Load spectrum data from the JSON fixture."""
    data = _FIXTURES[key]
    return {
        "spectrum": [tuple(p) for p in data["spectrum"]],
        "precursor_mz": data["precursor_mz"],
        "precursor_charge": data["precursor_charge"],
    }


class TestUSICase1(unittest.TestCase):
    """Known compound (scan:12430) with SMILES vs modified compound (scan:5538)."""

    KNOWN_SMILES = "CCC(C)(C)C(=O)O[C@H]1C[C@H](C=C2[C@H]1[C@H]([C@H](C=C2)C)CC[C@@H]3C[C@H](CC(=O)O3)O)C"

    @classmethod
    def setUpClass(cls):
        known_data = _load_fixture("case1_known")
        unknown_data = _load_fixture("case1_unknown")
        cls.known = Compound(
            spectrum=known_data["spectrum"],
            precursor_mz=known_data["precursor_mz"],
            precursor_charge=known_data["precursor_charge"],
            adduct="[M+H]1+",
            smiles=cls.KNOWN_SMILES,
            id="case1_known",
        )
        cls.unknown = Compound(
            spectrum=unknown_data["spectrum"],
            precursor_mz=unknown_data["precursor_mz"],
            precursor_charge=unknown_data["precursor_charge"],
            adduct="[M+H]1+",
            id="case1_unknown",
        )
        cls.mf = ModiFinder(
            knownCompound=cls.known,
            unknownCompound=cls.unknown,
            helpers=None,
            ppm_tolerance=40,
        )

    def test_compounds_loaded(self):
        self.assertIsNotNone(self.known.spectrum)
        self.assertIsNotNone(self.unknown.spectrum)
        self.assertIsNotNone(self.known.structure)

    def test_network_structure(self):
        self.assertEqual(len(self.mf.network.nodes), 2)
        self.assertEqual(len(self.mf.network.edges), 1)

    def test_generate_probabilities(self):
        scores = self.mf.generate_probabilities()
        self.assertEqual(len(scores), self.known.structure.GetNumAtoms())
        self.assertAlmostEqual(np.sum(scores), 1.0, places=5)
        self.assertTrue(np.all(scores >= 0))

    def test_modification_on_steroid_rings(self):
        """Modification site is on the steroid ring system (atoms 8-18, 29).
        The top predicted atom should be in that region."""
        scores = self.mf.generate_probabilities()
        predicted_site = int(np.argmax(scores))
        steroid_ring_atoms = list(range(8, 19)) + [29]
        self.assertIn(
            predicted_site,
            steroid_ring_atoms,
            f"Predicted site {predicted_site} is not on the steroid rings (atoms 8-18, 29)",
        )


class TestUSICase2(unittest.TestCase):
    """Known compound (scan:1850) with SMILES vs modified compound (scan:1728)."""

    KNOWN_SMILES = "CCCCCCCCCCCC=CC(=O)N(CCCCCNC(=O)CCC(=O)N(CCCCCNC(=O)C1COC(=N1)C2=CC=CC=C2O)O)O"

    @classmethod
    def setUpClass(cls):
        known_data = _load_fixture("case2_known")
        unknown_data = _load_fixture("case2_unknown")
        cls.known = Compound(
            spectrum=known_data["spectrum"],
            precursor_mz=known_data["precursor_mz"],
            precursor_charge=known_data["precursor_charge"],
            adduct="[M+H]1+",
            smiles=cls.KNOWN_SMILES,
            id="case2_known",
        )
        cls.unknown = Compound(
            spectrum=unknown_data["spectrum"],
            precursor_mz=unknown_data["precursor_mz"],
            precursor_charge=unknown_data["precursor_charge"],
            adduct="[M+H]1+",
            id="case2_unknown",
        )
        cls.mf = ModiFinder(
            knownCompound=cls.known,
            unknownCompound=cls.unknown,
            helpers=None,
            ppm_tolerance=40,
        )

    def test_compounds_loaded(self):
        self.assertIsNotNone(self.known.spectrum)
        self.assertIsNotNone(self.unknown.spectrum)
        self.assertIsNotNone(self.known.structure)

    def test_network_structure(self):
        self.assertEqual(len(self.mf.network.nodes), 2)
        self.assertEqual(len(self.mf.network.edges), 1)

    def test_generate_probabilities(self):
        scores = self.mf.generate_probabilities()
        self.assertEqual(len(scores), self.known.structure.GetNumAtoms())
        self.assertAlmostEqual(np.sum(scores), 1.0, places=5)
        self.assertTrue(np.all(scores >= 0))

    def test_top_prediction_on_tail(self):
        """The modification is on the alkyl tail (atoms 0-9).
        The top predicted atom should be in that region."""
        scores = self.mf.generate_probabilities()
        predicted_site = int(np.argmax(scores))
        tail_atoms = list(range(0, 10))
        self.assertIn(
            predicted_site,
            tail_atoms,
            f"Predicted site {predicted_site} is not on the alkyl tail (atoms 0-9)",
        )

    def test_get_edge_detail_no_side_effect(self):
        """Previously, calling get_edge_detail should not mutate the original edge matches.

        Bug description: EdgeDetail.copy() used a shallow copy of self.matches, so reverse_match()
        on the copy mutates the original Match objects.
        """
        uid = self.unknown.id
        kid = self.known.id

        if self.mf.network.has_edge(uid, kid):
            orig_edge = self.mf.network[uid][kid]["edgeDetail"]
        else:
            orig_edge = self.mf.network[kid][uid]["edgeDetail"]

        orig_first = [m.first_peak_mz for m in orig_edge.matches]
        orig_second = [m.second_peak_mz for m in orig_edge.matches]

        self.mf.get_edge_detail(kid, uid)

        after_first = [m.first_peak_mz for m in orig_edge.matches]
        after_second = [m.second_peak_mz for m in orig_edge.matches]

        self.assertEqual(orig_first, after_first, "first_peak_mz was mutated by get_edge_detail")
        self.assertEqual(orig_second, after_second, "second_peak_mz was mutated by get_edge_detail")


if __name__ == "__main__":
    unittest.main()
