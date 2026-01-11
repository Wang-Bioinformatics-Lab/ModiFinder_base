import unittest
from pathlib import Path

import numpy as np

from modifinder import Compound, ModiFinder, Spectrum
from modifinder.utilities.general_utils import read_mgf
from modifinder.utilities import mol_utils as mf_mu
from modifinder.utilities import visualizer as mf_vis


class TestPairMGFFlow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mgf_path = Path(__file__).resolve().parent / "test_files" / "pair.mgf"
        cls.mgf_df = read_mgf(str(cls.mgf_path))

    @staticmethod
    def _coerce_charge(charge):
        if isinstance(charge, (list, tuple)):
            return int(charge[0]) if charge else None
        if charge is None:
            return None
        return int(charge)

    @classmethod
    def _compound_from_row(cls, row):
        row = row.to_dict()
        charge = cls._coerce_charge(row.get("charge"))
        spectrum = row["spectrum"]
        data = {
            "mz": spectrum[0],
            "intensity": spectrum[1],
            "precursor_mz": row["precursor_mz"],
            "precursor_charge": charge,
            "smiles": row["smiles"],
            "name": row.get("name"),
            "spectrum_id": row.get("spectrumid"),
        }
        return Compound(data)

    def test_read_mgf(self):
        self.assertEqual(len(self.mgf_df), 2)
        self.assertIn("spectrum", self.mgf_df.columns)
        self.assertIn("smiles", self.mgf_df.columns)
        self.assertEqual(self.mgf_df.iloc[0]["spectrum"].shape[0], 2)

    def test_create_spectrum(self):
        row = self.mgf_df.iloc[0].to_dict()
        charge = self._coerce_charge(row.get("charge"))
        spectrum = Spectrum(
            mz=row["spectrum"][0],
            intensity=row["spectrum"][1],
            precursor_mz=row["precursor_mz"],
            precursor_charge=charge,
        )
        self.assertIsNotNone(spectrum.mz)
        self.assertEqual(len(spectrum.mz), 3)
        self.assertEqual(spectrum.precursor_charge, 1)

    def test_create_compound(self):
        compound = self._compound_from_row(self.mgf_df.iloc[0])
        self.assertIsNotNone(compound.structure)
        self.assertIsNotNone(compound.spectrum)
        symbols = {atom.GetSymbol() for atom in compound.structure.GetAtoms()}
        self.assertIn("N", symbols)

    def test_modifinder_solve_and_modification_site(self):
        known = self._compound_from_row(self.mgf_df.iloc[0])
        unknown = self._compound_from_row(self.mgf_df.iloc[1])
        mf = ModiFinder(known, unknown, ppm_tolerance=40)
        self.assertIsNone(mf.solve(unknown.id))
        mod_sites = mf_mu.get_modification_nodes(known.structure, unknown.structure, True)
        self.assertEqual(len(mod_sites), 1)
        n_indices = [atom.GetIdx() for atom in known.structure.GetAtoms() if atom.GetSymbol() == "N"]
        self.assertEqual(mod_sites[0], n_indices[0])

    def test_draw_spectrums(self):
        known = self._compound_from_row(self.mgf_df.iloc[0])
        unknown = self._compound_from_row(self.mgf_df.iloc[1])
        img_known = mf_vis.draw_spectrum(known.spectrum)
        img_unknown = mf_vis.draw_spectrum(unknown.spectrum)
        self.assertTrue(hasattr(img_known, "shape"))
        self.assertTrue(hasattr(img_unknown, "shape"))
        self.assertGreater(img_known.size, 0)
        self.assertGreater(img_unknown.size, 0)

    def test_draw_molecules(self):
        known = self._compound_from_row(self.mgf_df.iloc[0])
        unknown = self._compound_from_row(self.mgf_df.iloc[1])
        img_known = mf_vis.draw_molecule(known.structure)
        img_unknown = mf_vis.draw_molecule(unknown.structure)
        self.assertTrue(hasattr(img_known, "shape"))
        self.assertTrue(hasattr(img_unknown, "shape"))
        self.assertGreater(img_known.size, 0)
        self.assertGreater(img_unknown.size, 0)

    def test_draw_alignment(self):
        known = self._compound_from_row(self.mgf_df.iloc[0])
        unknown = self._compound_from_row(self.mgf_df.iloc[1])
        mf = ModiFinder(known, unknown, ppm_tolerance=40)
        img = mf.draw_alignment(known.id, unknown.id)
        self.assertTrue(hasattr(img, "shape"))
        self.assertGreater(img.size, 0)

    def test_draw_heatmap(self):
        known = self._compound_from_row(self.mgf_df.iloc[0])
        unknown = self._compound_from_row(self.mgf_df.iloc[1])
        mf = ModiFinder(known, unknown, ppm_tolerance=40)
        n_indices = [atom.GetIdx() for atom in known.structure.GetAtoms() if atom.GetSymbol() == "N"]
        scores = np.zeros(known.structure.GetNumAtoms())
        scores[n_indices[0]] = 1.0
        img = mf.draw_prediction(scores, known.id, show_legend=False)
        self.assertTrue(hasattr(img, "shape"))
        self.assertGreater(img.size, 0)


if __name__ == "__main__":
    unittest.main()
