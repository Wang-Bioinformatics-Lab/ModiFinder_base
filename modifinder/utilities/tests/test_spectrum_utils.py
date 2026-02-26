import unittest
from modifinder.utilities.spectrum_utils import aggregate_spectrums, refine_consensus
from modifinder import Spectrum

class TestIsShifted(unittest.TestCase):
    def __init__(self, methodName = "runTest"):
        super().__init__(methodName)
        self.positive_spectrums = []
        self.positive_spectrums.append(Spectrum(mz = [100, 200, 300], intensity = [1, 2, 3], adduct = "[M+H]+", precursor_mz = 1000))
        self.positive_spectrums.append(Spectrum(mz = [100, 199, 300], intensity = [1, 2, 3], adduct = "[M+H]+", precursor_mz = 1000))
        self.positive_spectrums.append(Spectrum(mz = [99, 200, 301], intensity = [1, 2, 3], adduct = "[M+H]+", precursor_mz = 1000))
        self.positive_spectrums.append(Spectrum(mz = [100, 200, 299], intensity = [1, 2, 3], adduct = "[M+H]+", precursor_mz = 1000))
        
        self.negative_spectrums = []
        h = 100 - 2 * 1.0078250321
        w = 200 - 2 * 1.0078250321
        t = 300 - 2 * 1.0078250321
        self.negative_spectrums.append(Spectrum(mz = [h, w, t], intensity = [1, 2, 3], adduct = "[M-H]-", precursor_mz = 1000))
        self.negative_spectrums.append(Spectrum(mz = [h, w - 1, t], intensity = [1, 2, 3], adduct = "[M-H]-", precursor_mz = 1000))
        self.negative_spectrums.append(Spectrum(mz = [h - 1, w, t + 1], intensity = [1, 2, 3], adduct = "[M-H]-", precursor_mz = 1000))
        self.negative_spectrums.append(Spectrum(mz = [h, w, t - 1], intensity = [1, 2, 3], adduct = "[M-H]-", precursor_mz = 1000))
        
    
    def refine_consensus(self):
        spectrums = self.positive_spectrums + self.negative_spectrums
        spectrum = spectrums[0]
        spectrums = spectrums[1:]
        result, appearances = refine_consensus(spectrum, spectrums, ppm_tolerance=1000,
                                               mz_tolerance=1000, consensus_majority_ratio = 0.5)
        self.assertEqual(result.mz, [100, 200, 300])
        self.assertEqual(result.intensity, [1, 2, 3])
        self.assertEqual(result.adduct, "[M+H]+")
        self.assertEqual(result.precursor_mz, 1000)
        

if __name__ == '__main__':
    unittest.main()