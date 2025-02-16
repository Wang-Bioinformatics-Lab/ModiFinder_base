
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors
from modifinder.utilities.network import get_np_classifier
from rdkit import Chem
from dataclasses import dataclass

@dataclass(frozen=True)
class StructureMeta:
    """
    Class to store metadata about a molecule.
    
    Attributes
    ----------
    num_aromatic_rings : int
        Number of aromatic rings in the molecule.
    num_atoms : int
        Number of atoms in the molecule.
    num_bonds : int
        Number of bonds in the molecule.
    num_rings : int
        Number of rings in the molecule.
    class_results: list
        List of classes from the np_classifier
    superclass: list
        List of superclasses from the np_classifier
    pathway: list
        List of pathways from the np_classifier
    isglycoside: bool
        True if the molecule is a glycoside.
    """
    
    num_aromatic_rings: int
    num_atoms: int
    num_bonds: int
    num_rings: int
    class_types: list
    superclasses: list
    pathways: list
    isglycoside: bool
    
    
    @classmethod
    def from_structure(cls, structure):
        """
        Create a StructureMeta object from an RDKit molecule.
        
        Parameters
        ----------
        structure : RDKit molecule
            RDKit molecule object.
        
        Returns
        -------
        StructureMeta
            StructureMeta object with metadata about the molecule.
        """
        
        num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(structure)
        num_atoms = structure.GetNumAtoms()
        num_bonds = structure.GetNumBonds()
        num_rings = rdMolDescriptors.CalcNumRings(structure)
        np_data = get_np_classifier(Chem.MolToSmiles(structure))
        class_types = np_data["class_types"]
        superclasses = np_data["superclasses"]
        pathways = np_data["pathways"]
        isglycoside = np_data["isglycoside"]
        return cls(num_aromatic_rings, num_atoms, num_bonds, num_rings, class_types, superclasses, pathways, isglycoside)