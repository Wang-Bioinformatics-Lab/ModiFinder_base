# Getting Started with ModiFinder

ModiFinder is your Swiss Army knife for mass spectrometry data analysis in Python. This guide will get you up and running quickly.

## Installation
ModiFinder requires Python 3.9 or above. 

This turorial assumes you have already setup the Python environment you want to work in and intend to install ``modifinder`` inside of it.  If you want
to create and work with Python virtual environments, please follow instructions
on `venv <https://docs.python.org/3/library/venv.html>`_ and `virtual
environments <http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_ or use `Conda/Mamba <https://github.com/mamba-org/mamba>`_ to create a new environment.

First, make sure you have the latest version of ``pip`` installed. `Pip documentation

<https://pip.pypa.io/en/stable/installing/>`_ 
### Using pip (Recommended)

```bash
pip install modifinder
```

### From Source

```bash
git clone https://github.com/your-repo/modifinder.git
cd modifinder
pip install -e .
```

### Verify Installation

```python
import modifinder
print(modifinder.__version__)
```

## Core Concepts

ModiFinder revolves around two main classes:

- **Spectrum**: Represents MS/MS spectral data
- **Compound**: Represents a chemical compound with its structure and spectrum

These classes can be created from various sources: GNPS identifiers, raw data arrays, MGF files, SMILES strings, and more.

## Your First Steps

### 1. Fetch and Visualize a Compound

```python
from modifinder import Compound
from modifinder.utilities import visualizer as viz
import matplotlib.pyplot as plt

# Get compound from GNPS
compound = Compound("CCMSLIB00010113829")

# Draw the molecule
img = viz.draw_molecule(compound.structure, label=compound.name)

plt.figure(figsize=(8, 8))
plt.imshow(img)
plt.axis('off')
plt.show()

print(f"Compound: {compound.name}")
print(f"Formula: {compound.formula if hasattr(compound, 'formula') else 'N/A'}")
print(f"Peaks: {len(compound.spectrum.mz)}")
```

### 2. Work with Spectra

```python
from modifinder import Spectrum

# Create from GNPS
spectrum = Spectrum("mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00010113829")

# Access data
print(f"Precursor m/z: {spectrum.precursor_mz}")
print(f"Number of peaks: {len(spectrum.mz)}")
print(f"Base peak m/z: {spectrum.mz[spectrum.intensity.argmax()]}")

# Normalize
spectrum.normalize()
print(f"Max intensity after normalization: {max(spectrum.intensity)}")
```

### 3. Process MGF Files

```python
from modifinder.utilities import general_utils as gu

# Read MGF file
df = gu.read_mgf("data.mgf")

print(f"Loaded {len(df)} spectra")
print(f"Columns: {df.columns.tolist()}")

# Filter by precursor mass
filtered = df[(df['precursor_mz'] >= 200) & (df['precursor_mz'] <= 500)]
print(f"Filtered to {len(filtered)} spectra")
```

### 4. Compare Structures

```python
from modifinder import Compound
from modifinder.utilities import visualizer as viz

# Two related compounds
compound1 = Compound("CCMSLIB00010113829")
compound2 = Compound("CCMSLIB00010125628")

# Visualize the difference
img = viz.draw_modifications(
    compound1.structure,
    compound2.structure,
    show_legend=True
)

plt.figure(figsize=(12, 6))
plt.imshow(img)
plt.axis('off')
plt.title("Structural Comparison")
plt.show()
```

### 5. Find Modification Sites

Now for ModiFinder's main purpose:

```python
from modifinder import ModiFinder, Compound

# Known and modified compounds
known = Compound("CCMSLIB00010113829")
modified = Compound("CCMSLIB00010125628")

# Run ModiFinder
mf = ModiFinder(known, modified, mz_tolerance=0.01, ppm_tolerance=40)
probabilities = mf.generate_probabilities()

# Visualize prediction
img = mf.draw_prediction(
    probabilities,
    known.id,
    show_legend=True,
    show_labels=True
)

plt.figure(figsize=(10, 10))
plt.imshow(img)
plt.axis('off')
plt.title("Predicted Modification Sites")
plt.show()
```

## Common Use Cases

### Extract Data from GNPS

```python
from modifinder.utilities import network

# Fetch data
data = network.get_data("CCMSLIB00010113829")

print(f"Name: {data['compound_name']}")
print(f"SMILES: {data['smiles']}")
print(f"Precursor: {data['precursor_mz']}")
```

### Batch Process Multiple Compounds

```python
from modifinder import Compound

accessions = ["CCMSLIB00010113829", "CCMSLIB00010125628", "CCMSLIB00010114304"]

compounds = []
for acc in accessions:
    try:
        compound = Compound(acc)
        compounds.append(compound)
        print(f"✓ {acc}: {compound.name}")
    except Exception as e:
        print(f"✗ {acc}: {e}")

print(f"\nSuccessfully loaded {len(compounds)} compounds")
```

### Create Publication-Quality Figures

```python
from modifinder import Compound
from modifinder.utilities import visualizer as viz

compound = Compound("CCMSLIB00010113829")

# High-resolution molecule image
img_mol = viz.draw_molecule(
    compound.structure,
    size=(1200, 1200),
    label=compound.name,
    label_font_size=48
)

# High-resolution spectrum
img_spec = viz.draw_spectrum(
    compound.spectrum,
    title=f"MS/MS Spectrum - m/z {compound.precursor_mz:.2f}"
)

# Save high-DPI
from PIL import Image
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
axes[0].imshow(img_mol)
axes[0].axis('off')
axes[1].imshow(img_spec)
axes[1].axis('off')
plt.savefig("publication_figure.png", dpi=300, bbox_inches='tight')
```

## What's Next?

Explore the detailed tutorials:

- **[Working with USIs](tutorials/working_with_usi.md)**: Fetch and process data from GNPS
- **[Working with MGF Files](tutorials/working_with_mgf.md)**: Read, filter, and batch process MGF files
- **[Visualization](tutorials/visualization.md)**: Create beautiful figures of molecules and spectra
- **[ModiFinder Basics](tutorials/basics.md)**: Deep dive into modification site prediction
- **[Customization](tutorials/customization.md)**: Advanced ModiFinder configuration

## Need Help?

- Check the [API Reference](modifinder/index.rst) for detailed function documentation
- See example notebooks in the `docs/source/tutorials/` directory
- Report issues on GitHub

## Tips for Success

1. **Start Simple**: Begin with GNPS identifiers before working with custom data
2. **Validate Data**: Always check that objects were created successfully
3. **Visualize Often**: Use the drawing tools to understand your data
4. **Process in Batches**: For large datasets, process data in manageable chunks
5. **Preserve Metadata**: Keep track of IDs and other metadata throughout your workflow

Happy analyzing! 🚀
