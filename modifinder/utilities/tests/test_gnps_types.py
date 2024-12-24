import unittest

from modifinder.utilities.gnps_types import (
    convert_to_SpectrumTuple,
    convert_to_SpectrumTuple_seprated,
    convert_to_universal_key,
    parse_data_to_universal,
    Convert_SpectrumTuple_to_peaks,
    SpectrumTuple
)

class TestGNPSTypes(unittest.TestCase):

    def test_convert_to_SpectrumTuple(self):
        peaks = [(100.0, 200.0), (150.0, 300.0)]
        precursor_mz = 500.0
        precursor_charge = 2
        expected = SpectrumTuple(precursor_mz=500.0, precursor_charge=2, mz=[100.0, 150.0], intensity=[200.0, 300.0])
        result = convert_to_SpectrumTuple(peaks, precursor_mz, precursor_charge)
        self.assertEqual(result, expected)

        with self.assertRaises(ValueError):
            convert_to_SpectrumTuple([100.0, 200.0], precursor_mz, precursor_charge)

        self.assertIsNone(convert_to_SpectrumTuple([], precursor_mz, precursor_charge))

    def test_convert_to_SpectrumTuple_seprated(self):
        mz = [100.0, 150.0]
        intensity = [200.0, 300.0]
        precursor_mz = 500.0
        precursor_charge = 2
        expected = SpectrumTuple(precursor_mz=500.0, precursor_charge=2, mz=[100.0, 150.0], intensity=[200.0, 300.0])
        result = convert_to_SpectrumTuple_seprated(mz, intensity, precursor_mz, precursor_charge)
        self.assertEqual(result, expected)

        with self.assertRaises(ValueError):
            convert_to_SpectrumTuple_seprated([100.0], [200.0, 300.0], precursor_mz, precursor_charge)

        self.assertIsNone(convert_to_SpectrumTuple_seprated([], [], precursor_mz, precursor_charge))

    def test_convert_to_universal_key(self):
        self.assertEqual(convert_to_universal_key("precursor_mz"), "Precursor_MZ")
        self.assertEqual(convert_to_universal_key("smiles"), "Smiles")
        self.assertEqual(convert_to_universal_key("SMILES"), "Smiles")
        self.assertEqual(convert_to_universal_key("charge"), "Charge")
        self.assertEqual(convert_to_universal_key("adduct"), "Adduct")
        self.assertEqual(convert_to_universal_key("unknown_key"), "unknown_key")

    def test_parse_data_to_universal(self):
        data = {
            "peaks_json": '[{"mz": 100.0, "intensity": 200.0}, {"mz": 150.0, "intensity": 300.0}]',
            "precursor_mz": "500.0",
            "Charge": "2"
        }
        expected = {
            "peaks": [{"mz": 100.0, "intensity": 200.0}, {"mz": 150.0, "intensity": 300.0}],
            "Precursor_MZ": 500.0,
            "Charge": 2
        }
        result = parse_data_to_universal(data)
        self.assertEqual(result, expected)

    def test_Convert_SpectrumTuple_to_peaks(self):
        spectrum = SpectrumTuple(precursor_mz=500.0, precursor_charge=2, mz=[100.0, 150.0], intensity=[200.0, 300.0])
        expected = [(100.0, 200.0), (150.0, 300.0)]
        result = Convert_SpectrumTuple_to_peaks(spectrum)
        self.assertEqual(result, expected)

if __name__ == '__main__':
    unittest.main()