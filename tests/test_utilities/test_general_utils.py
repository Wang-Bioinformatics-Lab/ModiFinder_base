import unittest
from modifinder.utilities.general_utils import is_shifted, parse_data_to_universal

class TestIsShifted(unittest.TestCase):
    def test_is_shifted_with_ppm(self):
        self.assertTrue(is_shifted(100.0, 100.1, ppm_tolerance=500))
        self.assertFalse(is_shifted(100.0, 100.0001, ppm_tolerance=500))

    def test_is_shifted_with_mz_tol(self):
        self.assertTrue(is_shifted(100.0, 100.2, mz_tolerance=0.1))
        self.assertFalse(is_shifted(100.0, 100.05, mz_tolerance=0.1))

    def test_is_shifted_with_both_ppm_and_mz_tol(self):
        self.assertTrue(is_shifted(100.0, 100.2, ppm_tolerance=500, mz_tolerance=0.1))
        self.assertTrue(is_shifted(100.0, 100.0002, ppm_tolerance=1, mz_tolerance=0.0001))
        self.assertFalse(is_shifted(100.0, 100.00005, ppm_tolerance=1, mz_tolerance=0.0001))

    def test_is_shifted_without_ppm_or_mz_tol(self):
        with self.assertRaises(ValueError):
            is_shifted(100.0, 100.1)
    
    def test_parse_data_to_universal(self):
        data = {
            "peaks_json": '[{"mz": 100.0, "intensity": 200.0}, {"mz": 150.0, "intensity": 300.0}]',
            "precursor_mz": "500.0",
            "Charge": "2"
        }
        expected = {
            "peaks": [{"mz": 100.0, "intensity": 200.0}, {"mz": 150.0, "intensity": 300.0}],
            "precursor_mz": 500.0,
            "precursor_charge": 2
        }
        result = parse_data_to_universal(data)
        self.assertEqual(result, expected)

if __name__ == '__main__':
    unittest.main()
