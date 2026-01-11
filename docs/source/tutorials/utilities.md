# ModiFinder Utilities: Essential Tools for MS Data Processing

ModiFinder includes powerful utility functions for common mass spectrometry data tasks. This tutorial covers the most useful tools for spectrum processing, network data retrieval, and molecule manipulation.

## Working with Network Data

ModiFinder can fetch spectral data directly from GNPS using USI or accession numbers.

### Fetching Spectrum Data

```python
from modifinder.utilities.network import get_data

# Using GNPS accession
data = get_data("CCMSLIB00005436077")
print(f"Precursor m/z: {data['precursor_mz']}")
print(f"Charge: {data['precursor_charge']}")
print(f"Number of peaks: {len(data['mz'])}")

# Using USI
usi = "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00005436077"
data_from_usi = get_data(usi)
```

### Getting Matched Peaks

Compare two spectra and get their matched peaks using the modified cosine algorithm:

```python
from modifinder.utilities.network import get_matched_peaks

# Compare two spectra
identifier1 = "CCMSLIB00005436077"
identifier2 = "CCMSLIB00005436078"

matches = get_matched_peaks(identifier1, identifier2)
print(f"Cosine score: {matches.get('cosine', 'N/A')}")
```

### Natural Product Classification

Get chemical classification information for a molecule:

```python
from modifinder.utilities.network import get_np_classifier

smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # caffeine
classification = get_np_classifier(smiles)

print(f"Pathways: {classification['pathways']}")
print(f"Superclasses: {classification['superclasses']}")
print(f"Is glycoside: {classification['isglycoside']}")
```

## Consensus Spectra Creation

When you have multiple spectra of the same compound, create a consensus spectrum to reduce noise and improve quality.

### Aggregating Multiple Spectra

```python
from modifinder.utilities.spectrum_utils import aggregate_spectrums
from modifinder import Spectrum

# Load multiple spectra of the same compound
spectrums = [
    Spectrum("mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00005436077"),
    Spectrum("mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00005436078"),
    # Add more spectra...
]

# Create consensus spectrum
consensus = aggregate_spectrums(
    spectrums,
    ppm_tolerance=10,
    mz_tolerance=0.1,
    consensus_majority_ratio=0.5  # Peak must appear in 50% of spectra
)

print(f"Original peaks: {len(spectrums[0].mz)}")
print(f"Consensus peaks: {len(consensus.mz)}")
```

### Refining a Spectrum

Remove noisy peaks from a spectrum by comparing against reference spectra:

```python
from modifinder.utilities.spectrum_utils import refine_consensus

# Refine a spectrum based on consensus from other spectra
refined, (new_appearances, old_appearances) = refine_consensus(
    spectrum=spectrums[0],
    spectrums=spectrums[1:],
    ppm_tolerance=10,
    mz_tolerance=0.1,
    consensus_majority_ratio=0.3  # Keep peaks appearing in 30% of references
)

print(f"Peaks before refinement: {len(spectrums[0].mz)}")
print(f"Peaks after refinement: {len(refined.mz)}")
```

## Working with MGF Files

Read and write MGF (Mascot Generic Format) files easily:

### Reading MGF Files

```python
from modifinder.utilities.general_utils import read_mgf

# Read MGF file into DataFrame
df = read_mgf("path/to/spectra.mgf")

# Access spectrum data
for idx, row in df.iterrows():
    mz_array = row['spectrum'][0]
    intensity_array = row['spectrum'][1]
    precursor = row['precursor_mz']
    print(f"Spectrum {idx}: {len(mz_array)} peaks, precursor={precursor}")
```

### Writing MGF Files

```python
from modifinder.utilities.general_utils import write_mgf
import pandas as pd
import numpy as np

# Create or modify DataFrame
msms_df = pd.DataFrame({
    'precursor_mz': [500.123, 600.456],
    'charge': [2, 1],
    'spectrum': [
        np.array([[100, 200], [50, 100]]),  # [mz_array, intensity_array]
        np.array([[150, 250], [75, 125]])
    ]
})

write_mgf(msms_df, "output.mgf")
```

## Molecular Formula Utilities

### Getting Formula Mass

```python
from modifinder.utilities.general_utils import get_formula_mass

mass = get_formula_mass("C8H10N4O2")  # caffeine
print(f"Molecular mass: {mass:.4f}")
```

### Getting Adduct Mass

```python
from modifinder.utilities.general_utils import get_adduct_mass

# Standard adducts
adduct_mass = get_adduct_mass("[M+H]+")
print(f"[M+H]+ mass shift: {adduct_mass:.4f}")

adduct_mass = get_adduct_mass("[M+Na]+")
print(f"[M+Na]+ mass shift: {adduct_mass:.4f}")

# More complex adducts
adduct_mass = get_adduct_mass("[2M+H]+")
print(f"[2M+H]+ mass shift: {adduct_mass:.4f}")
```

### Parsing Molecular Formulas

```python
from modifinder.utilities.general_utils import parse_molecular_formula

formula_dict = parse_molecular_formula("C8H10N4O2")
print(formula_dict)  # {'C': 8, 'H': 10, 'N': 4, 'O': 2}
```

### Checking Substructures

```python
from modifinder.utilities.general_utils import is_submolecule

# Check if one formula is a substructure of another
is_sub = is_submolecule("C6H12O6", "C12H22O11", ignore_H=False)
print(f"Is glucose part of sucrose? {is_sub}")
```

## Molecule Utilities

### Calculating Edit Distance

Measure structural similarity between two molecules:

```python
from modifinder.utilities.mol_utils import get_edit_distance
from rdkit import Chem

mol1 = Chem.MolFromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")  # caffeine
mol2 = Chem.MolFromSmiles("CN1C=NC2=C1C(=O)NC(=O)N2C")     # theophylline

distance = get_edit_distance(mol1, mol2)
print(f"Edit distance: {distance}")
```

### Finding Modification Sites

Identify which atoms differ between two molecules:

```python
from modifinder.utilities.mol_utils import get_modification_nodes
from rdkit import Chem

mol1 = Chem.MolFromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
mol2 = Chem.MolFromSmiles("CN1C=NC2=C1C(=O)NC(=O)N2C")

# Get modification sites in mol1
mod_sites_mol1 = get_modification_nodes(mol1, mol2, in_mol1=True)
print(f"Modified atoms in mol1: {mod_sites_mol1}")

# Get modification sites in mol2
mod_sites_mol2 = get_modification_nodes(mol1, mol2, in_mol1=False)
print(f"Modified atoms in mol2: {mod_sites_mol2}")
```

### Getting Transition Details

Get comprehensive information about the structural changes between two molecules:

```python
from modifinder.utilities.mol_utils import get_transition

transition = get_transition(mol1, mol2)

print(f"Common atoms: {transition['common_atoms']}")
print(f"Removed atoms: {transition['removed_atoms']}")
print(f"Added atoms: {transition['added_atoms']}")
print(f"Modified bonds (added): {len(transition['modified_added_edges_inside'])}")
print(f"Modified bonds (removed): {len(transition['modified_removed_edges_inside'])}")
```

## Tolerance and Comparison Utilities

### Checking if Values are Shifted

Determine if two m/z values differ beyond tolerance:

```python
from modifinder.utilities.general_utils import is_shifted

# Using ppm tolerance
shifted = is_shifted(500.0, 500.01, ppm_tolerance=10)
print(f"Values shifted (10 ppm)? {shifted}")

# Using absolute m/z tolerance
shifted = is_shifted(500.0, 500.15, mz_tolerance=0.1)
print(f"Values shifted (0.1 Da)? {shifted}")

# Using both (returns True if either condition is met)
shifted = is_shifted(500.0, 500.15, ppm_tolerance=10, mz_tolerance=0.1)
print(f"Values shifted (either condition)? {shifted}")
```

## Summary

These utilities make ModiFinder a comprehensive toolkit for MS data analysis:

- **Network utilities**: Fetch data from GNPS, get spectral matches, classify natural products
- **Spectrum utilities**: Create consensus spectra, refine spectra, handle adducts
- **File I/O**: Read and write MGF files efficiently
- **Formula utilities**: Calculate masses, parse formulas, check substructures
- **Molecule utilities**: Calculate edit distances, find modification sites, analyze transitions
- **Comparison utilities**: Handle m/z tolerance checking with ppm or absolute values

These tools work seamlessly with ModiFinder's core functionality and can be used independently for various MS data processing tasks.
