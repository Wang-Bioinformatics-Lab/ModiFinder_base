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

    def test_tail_high_probability(self):
        """The modification is on the alkyl tail (atoms 0-9).
        Tail atoms should have high probability and non-tail atoms should have low probability."""
        scores = self.mf.generate_probabilities()
        tail_atoms = list(range(0, 10))
        non_tail_atoms = list(range(10, len(scores)))

        tail_sum = sum(scores[i] for i in tail_atoms)
        tail_mean = np.mean([scores[i] for i in tail_atoms])
        non_tail_mean = np.mean([scores[i] for i in non_tail_atoms])

        # Tail should hold the majority of probability mass
        self.assertGreater(tail_sum, 0.5,
            f"Tail probability sum {tail_sum:.4f} should be > 0.5")

        # Tail mean should be significantly higher than non-tail mean
        self.assertGreater(tail_mean, non_tail_mean * 3,
            f"Tail mean {tail_mean:.4f} should be much higher than non-tail mean {non_tail_mean:.4f}")

        # Every non-tail atom should have low probability
        for i in non_tail_atoms:
            self.assertLess(scores[i], 0.03,
                f"Non-tail atom {i} has unexpectedly high probability {scores[i]:.4f}")

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

    def test_get_result_returns_valid_dicts(self):
        """get_result() should return (peaksObj, fragmentsObj) with correct keys and types."""
        peaksObj, fragmentsObj = self.mf.get_result()

        # peaksObj should be a dict with the expected keys
        self.assertIsInstance(peaksObj, dict)
        expected_peaks_keys = [
            "main_compound_peaks", "mod_compound_peaks", "matched_peaks",
            "main_precursor_mz", "mod_precursor_mz",
        ]
        for key in expected_peaks_keys:
            self.assertIn(key, peaksObj, f"peaksObj missing key: {key}")

        # fragmentsObj should be a dict with the expected keys
        self.assertIsInstance(fragmentsObj, dict)
        expected_frags_keys = ["frags_map", "structure", "peaks", "Precursor_MZ"]
        for key in expected_frags_keys:
            self.assertIn(key, fragmentsObj, f"fragmentsObj missing key: {key}")

    def test_get_result_peaks_content(self):
        """get_result() peaks lists should be non-empty tuples of length 2."""
        peaksObj, _ = self.mf.get_result()

        self.assertGreater(len(peaksObj["main_compound_peaks"]), 0)
        self.assertGreater(len(peaksObj["mod_compound_peaks"]), 0)

        # Each peak should be a (mz, intensity) pair
        for peak in peaksObj["main_compound_peaks"]:
            self.assertEqual(len(peak), 2)
        for peak in peaksObj["mod_compound_peaks"]:
            self.assertEqual(len(peak), 2)

    def test_get_result_precursor_mz(self):
        """get_result() precursor_mz values should match the compounds."""
        peaksObj, fragmentsObj = self.mf.get_result()

        self.assertAlmostEqual(peaksObj["main_precursor_mz"],
                               self.known.spectrum.precursor_mz)
        self.assertAlmostEqual(peaksObj["mod_precursor_mz"],
                               self.unknown.spectrum.precursor_mz)
        self.assertAlmostEqual(fragmentsObj["Precursor_MZ"],
                               self.known.spectrum.precursor_mz)

    def test_get_result_fragments_structure(self):
        """fragmentsObj structure and peaks should reference the known compound."""
        from rdkit import Chem
        peaksObj, fragmentsObj = self.mf.get_result()

        # Compare by SMILES since RDKit Mol objects may be copied internally
        self.assertEqual(
            Chem.MolToSmiles(fragmentsObj["structure"]),
            Chem.MolToSmiles(self.known.structure),
        )
        self.assertEqual(fragmentsObj["peaks"], peaksObj["main_compound_peaks"])

    def test_get_result_matched_peaks(self):
        """get_result() matched_peaks should contain expected pairs (unknown_mz, known_mz)."""
        peaksObj, _ = self.mf.get_result()
        matched = peaksObj["matched_peaks"]

        self.assertIsInstance(matched, list)
        self.assertEqual(len(matched), 9)

        expected_matched_peaks = [
            (409304992, 313139007),
            (308160003, 308160003),
            (209190002, 113023002),
            (310273986, 214106994),
            (102091003, 102091003),
            (84081001, 84080001),
            (390165008, 390165985),
            (201123001, 201123001),
            (165102005, 165102005),
        ]

        self.assertEqual(matched, expected_matched_peaks,
                         "matched_peaks do not match expected values")


if __name__ == "__main__":
    unittest.main()
