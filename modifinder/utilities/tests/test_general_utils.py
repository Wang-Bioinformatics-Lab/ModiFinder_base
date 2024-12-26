import unittest
from modifinder.utilities.general_utils import is_shifted

class TestIsShifted(unittest.TestCase):
    def test_is_shifted_with_ppm(self):
        self.assertTrue(is_shifted(100.0, 100.1, ppm=500))
        self.assertFalse(is_shifted(100.0, 100.0001, ppm=500))

    def test_is_shifted_with_mz_tol(self):
        self.assertTrue(is_shifted(100.0, 100.2, mz_tol=0.1))
        self.assertFalse(is_shifted(100.0, 100.05, mz_tol=0.1))

    def test_is_shifted_with_both_ppm_and_mz_tol(self):
        self.assertTrue(is_shifted(100.0, 100.2, ppm=500, mz_tol=0.1))
        self.assertTrue(is_shifted(100.0, 100.0002, ppm=1, mz_tol=0.0001))
        self.assertFalse(is_shifted(100.0, 100.00005, ppm=1, mz_tol=0.0001))

    def test_is_shifted_without_ppm_or_mz_tol(self):
        with self.assertRaises(ValueError):
            is_shifted(100.0, 100.1)

if __name__ == '__main__':
    unittest.main()