"""Test cases using USI-based spectra from GNPS2.

Case 1: Known compound (scan:12430) with SMILES + Modified compound (scan:5538).
Case 2: Known compound (scan:1850) with SMILES + Modified compound (scan:1728).
         This is the typical dashboard use case.
"""

import unittest
import requests
import numpy as np

from modifinder import Compound, ModiFinder


def _fetch_usi(usi: str) -> dict:
    """Fetch spectrum data from a USI via the GNPS2 API."""
    url = "https://metabolomics-usi.gnps2.org/json/"
    r = requests.get(url, params={"usi1": usi}, timeout=30)
    r.raise_for_status()
    data = r.json()
    spectrum = [tuple(p) for p in data["peaks"]]
    return {
        "spectrum": spectrum,
        "precursor_mz": float(data["precursor_mz"]),
        "precursor_charge": int(data.get("precursor_charge", 1) or 1),
    }


class TestUSICase1(unittest.TestCase):
    """Known compound (scan:12430) with SMILES vs modified compound (scan:5538)."""

    KNOWN_USI = "mzspec:GNPS2:TASK-1b5b94b4191d4223a5f57afb2aaaf0b0-nf_output/clustering/spectra_reformatted.mgf:scan:12430"
    UNKNOWN_USI = "mzspec:GNPS2:TASK-1b5b94b4191d4223a5f57afb2aaaf0b0-nf_output/clustering/spectra_reformatted.mgf:scan:5538"
    KNOWN_SMILES = "CCC(C)(C)C(=O)O[C@H]1C[C@H](C=C2[C@H]1[C@H]([C@H](C=C2)C)CC[C@@H]3C[C@H](CC(=O)O3)O)C"

    @classmethod
    def setUpClass(cls):
        known_data = _fetch_usi(cls.KNOWN_USI)
        unknown_data = _fetch_usi(cls.UNKNOWN_USI)
        cls.known = Compound(
            spectrum=known_data["spectrum"],
            precursor_mz=known_data["precursor_mz"],
            precursor_charge=known_data["precursor_charge"],
            adduct="[M+H]1+",
            smiles=cls.KNOWN_SMILES,
            id=cls.KNOWN_USI,
        )
        cls.unknown = Compound(
            spectrum=unknown_data["spectrum"],
            precursor_mz=unknown_data["precursor_mz"],
            precursor_charge=unknown_data["precursor_charge"],
            adduct="[M+H]1+",
            id=cls.UNKNOWN_USI,
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

    KNOWN_USI = "mzspec:GNPS2:TASK-f121ff43e4e2424ca719d85de44fbf19-nf_output/clustering/specs_ms.mgf:scan:1850"
    UNKNOWN_USI = "mzspec:GNPS2:TASK-f121ff43e4e2424ca719d85de44fbf19-nf_output/clustering/specs_ms.mgf:scan:1728"
    KNOWN_SMILES = "CCCCCCCCCCCC=CC(=O)N(CCCCCNC(=O)CCC(=O)N(CCCCCNC(=O)C1COC(=N1)C2=CC=CC=C2O)O)O"

    @classmethod
    def setUpClass(cls):
        known_data = _fetch_usi(cls.KNOWN_USI)
        unknown_data = _fetch_usi(cls.UNKNOWN_USI)
        cls.known = Compound(
            spectrum=known_data["spectrum"],
            precursor_mz=known_data["precursor_mz"],
            precursor_charge=known_data["precursor_charge"],
            adduct="[M+H]1+",
            smiles=cls.KNOWN_SMILES,
            id=cls.KNOWN_USI,
        )
        cls.unknown = Compound(
            spectrum=unknown_data["spectrum"],
            precursor_mz=unknown_data["precursor_mz"],
            precursor_charge=unknown_data["precursor_charge"],
            adduct="[M+H]1+",
            id=cls.UNKNOWN_USI,
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

    @unittest.expectedFailure
    def test_get_edge_detail_no_side_effect(self):
        """Calling get_edge_detail should not mutate the original edge matches.

        Known bug: EdgeDetail.copy() uses shallow list copy, so reverse_match()
        on the copy mutates the original Match objects.
        """
        uid = self.unknown.id
        kid = self.known.id

        if self.mf.network.has_edge(uid, kid):
            orig_edge = self.mf.network[uid][kid]["edgedetail"]
        else:
            orig_edge = self.mf.network[kid][uid]["edgedetail"]

        orig_first = [m.first_peak_mz for m in orig_edge.matches]
        orig_second = [m.second_peak_mz for m in orig_edge.matches]

        self.mf.get_edge_detail(kid, uid)

        after_first = [m.first_peak_mz for m in orig_edge.matches]
        after_second = [m.second_peak_mz for m in orig_edge.matches]

        self.assertEqual(orig_first, after_first, "first_peak_mz was mutated by get_edge_detail")
        self.assertEqual(orig_second, after_second, "second_peak_mz was mutated by get_edge_detail")


if __name__ == "__main__":
    unittest.main()
