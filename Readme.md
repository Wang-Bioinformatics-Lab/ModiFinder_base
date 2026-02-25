Welcome to ModiFinder's documentation!
======================================

**ModiFinder** is a tool for site localization of structural modifications using MS/MS data.


_This project is under active development._

The documentation is available at: [https://wang-bioinformatics-lab.github.io/ModiFinder_base/](https://wang-bioinformatics-lab.github.io/ModiFinder_base/)


Citing
------

    ModiFinder: Tandem Mass Spectral Alignment Enables Structural Modification Site Localization

    Mohammad Reza Zare Shahneh, Michael Strobel, Giovanni Andrea Vitale, Christian Geibel, Yasin El Abiead, Neha Garg, Berenike Wagner, Karl Forchhammer, Allegra Aron, Vanessa V Phelan, Daniel Petras, and Mingxun Wang

    Journal of the American Society for Mass Spectrometry 2024 35 (11), 2564-2578

    DOI: 10.1021/jasms.4c00061


License
-------
   
    Academic Software License: © 2024 UCR (“Institution”). Academic or nonprofit researchers are permitted to use this Software (as defined below) subject to Paragraphs 1-4:

    1. Institution hereby grants to you free of charge, so long as you are an academic or nonprofit researcher, a nonexclusive license under Institution’s copyright ownership interest in this software and any derivative works made by you thereof (collectively, the “Software”) to use, copy, and make derivative works of the Software solely for educational or academic research purposes, and to distribute such Software free of charge to other academic or nonprofit researchers for their educational or academic research purposes, in all cases subject to the terms of this Academic Software License. Except as granted herein, all rights are reserved by Institution, including the right to pursue patent protection of the Software.
    2. Any distribution of copies of this Software -- including any derivative works made by you thereof -- must include a copy (including the copyright notice above), and be made subject to the terms, of this Academic Software License; failure by you to adhere to the requirements in Paragraphs 1 and 2 will result in immediate termination of the license granted to you pursuant to this Academic Software License effective as of the date you first used the Software.
    3. IN NO EVENT WILL INSTITUTION BE LIABLE TO ANY ENTITY OR PERSON FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF INSTITUTION HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. INSTITUTION SPECIFICALLY DISCLAIMS ANY AND ALL WARRANTIES, EXPRESS AND IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE IS PROVIDED “AS IS.” INSTITUTION HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS OF THIS SOFTWARE.
    4. Any academic or scholarly publication arising from the use of this Software or any derivative works thereof will include the following acknowledgment:  The Software used in this research was created by [INSERT AUTHOR NAMES] of UC Riverside. © 2024 UCR.

    Commercial entities: please contact mingxun.wang@cs.ucr.edu or tp@ucr.edu for licensing opportunities.

Core API (Working)
--------
ModiFinder requires two spectrum objects:
```
main_compound = Compound(
    spectrum=s1_peaks,                       # Formatted as [[mz, int], ...]
    precursor_mz=s1_prec_mz,                 # Float
    precursor_charge=s1_charge,              # Int
    adduct=s1_adduct,                        # Str
    smiles=s1_smiles                         # Str
)
mod_compound = Compound(
    spectrum=s2_peaks,
    precursor_mz=s2_prec_mz,
    precursor_charge=s2_charge,
    adduct=s2_adduct,
    smiles=s2_smiles                         # Optional
)
```
Pass these to a ModiFinder Object
```
helper_compounds = None
siteLocator = ModiFinder(
    main_compound,
    mod_compound,
    helpers=helper_compounds,
    **args)
```
Collect dictionary results:
```
peaksObj, fragmentsObj = siteLocator.get_result()
```
Dictionary format:
```
peaksObj = {
            "main_compound_peaks": main_compound_peaks,
            "mod_compound_peaks": mod_compound_peaks,
            "matched_peaks": matched_peaks,
            "main_precursor_mz": knownCompound.spectrum.precursor_mz,
            "mod_precursor_mz": unknownCompound.spectrum.precursor_mz,
        }
```
```
fragmentsObj = {
    "frags_map": knownCompound.spectrum.peak_fragment_dict,
    "structure": knownCompound.structure,
    "peaks": main_compound_peaks,
    "Precursor_MZ": knownCompound.spectrum.precursor_mz,
}
```

Useful Utility Functions
------------------------
ModiFinder includes several useful utility functions for mass spectrometry data analysis and visualization, exposed under `modifinder.utilities`.

**Reading MGF Files**
You can easily read MGF files into a pandas DataFrame using `read_mgf`.

```python
from modifinder.utilities import read_mgf

# Read MGF file
df = read_mgf("path/to/your/spectrum.mgf")
print(df.head())
```

**Visualization**
The `vis` module provides powerful visualization tools.

```python
from modifinder.utilities import vis
import matplotlib.pyplot as plt

# Draw a molecule
img = vis.draw_molecule("C1=CC=C(C=C1)O", output_type="png")
plt.imshow(img)
plt.show()

# Draw a spectrum
# spectrum_data is a list of (mz, intensity) tuples
spectrum_data = df.iloc[0]['spectrum'].T.tolist()
vis.draw_spectrum(spectrum_data)
plt.show()
```

Developer Guide
---------------

### Running Tests
To run the automated tests, ensure you have `pytest` installed. Then run:

```bash
pytest
```
or specifically for utilities:
```bash
pytest modifinder/utilities/tests/
```

### Release Process
ModiFinder uses GitHub Actions for automated releases and documentation deployment.

**Creating a New Release (PyPI)**
1. Create and push a new tag starting with `v` (e.g., `v1.2.3`).
   ```bash
   git tag v1.2.3
   git push origin v1.2.3
   ```
   The version is automatically derived from the tag using `setuptools_scm`.
2. This triggers the `pypi_publish.yml` workflow, which:
   - Builds the package.
   - Publishes it to PyPI (and TestPyPI).
   - Creates a GitHub Release.

**Updating Documentation**
Documentation is automatically rebuilt and deployed to GitHub Pages whenever changes are pushed to the `main` branch (via `documentation.yml`).
