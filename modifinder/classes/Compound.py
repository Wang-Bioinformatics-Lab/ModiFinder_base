

import modifinder.utilities.general_utils as general_utils
from modifinder.utilities.mol_utils import _get_molecule
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors
from rdkit import Chem
from typing import Dict
import numpy as np
import math


# import modifinder as mf
from modifinder.classes.Spectrum import Spectrum
from modifinder import convert as convert


class Compound:
    """ A class to represent a compound 

    The compound always has spectrum data.
    If the compound structure is known, it also has a structure property.
    
    Attributes
    ----------
    Known Compound Attributes:
    
        id (str) : The id of the compound
            
        spectrum (Spectrum) : Spectrum Tuple representing mz, intensity, Precursor mass and Charge for mass spectrumetry data
    
    Known Compound Attributes:

        structure (Chem.Mol): The structure of the compound
        
        distances (dict): A dictionary of distances between every pair of atoms in the compound
        
        peak_fragments_map (dict): A dictionary mapping peaks to fragments

    Other Attributes:
    
        adduct_mass (float): The adduct mass, derived from the Adduct
        
        is_known (bool): A boolean indicating whether the compound is known, derived from the structure
        
        usi (str): The USI of the compound
        
        name (str): The name of the compound
        
        accession (str): The accession of the compound
        
        library_membership (str): The GNPS library membership of the compound
        
        exact_mass (float): The exact mass of the compound
    
    Examples
    --------
    Create a compound by providing the necessary information:

    >>> compound = Compound(id="CCMSLIB00005435812", peaks=[[110.066925,38499.089844],[138.060638,412152.093750],[195.079575,6894530.000000],[195.200180,480874.812500],[196.082092,43027.628906]], precursor_mz=195.087, precursor_charge=1, adduct="[M+H]+", smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

    Alternatively, you can create a compound by providing a dictionary of data:

    >>> data = {
    ...     "id": "CCMSLIB00005435812",
    ...     "peaks": [[110.066925,38499.089844],[138.060638,412152.093750],[195.079575,6894530.000000],[195.200180,480874.812500],[196.082092,43027.628906]],
    ...     "precursor_mz": 195.087,
    ...     "charge": 1,
    ...     "adduct": "[M+H]+",
    ...     "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    ... }
    >>> compound = Compound(data)

    You can also create a compound by providing a usi:

    >>> usi = "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00005435812"
    >>> compound = Compound(usi)

    or with an accession:

    >>> accession = "CCMSLIB00005435812"
    >>> compound = Compound(accession)
    
    """

    # def __init__(self, data = None, structure = None, id: str = None, spectrum: Spectrum = None,
    #             is_known: bool = None, name: str = None, peak_fragments_map: dict = None, distances: dict = None, **kwargs):
    def __init__(self, incoming_data=None, **kwargs):
        """Initialize the Compound class

        The compound class can be initialized in three different ways:
        1. By providing a dictionary of data to the *data* parameter that contains all the necessary information
        2. By providing a usi to the *data* parameter to retrieve the necessary information from GNPS
        3. By providing the necessary information as parameter
        
        If both the *data* and the parameters are provided, the parameters will override the data.
        
        If of the incoming data and kwargs are provided, an empty instance of the class will be created.

        Parameters
        ----------
        incoming_data : input data (optional, default: None)
            The data to initialize the class with, can be a dictionary of data, a usi, or a compound object. If not provided,
            the class will be initialized with the provided keyword arguments. If provided, the keyword arguments will still
            override the data.

        kwargs : keyword arguments (optional, default: No attributes)
            Attributes to initialize the class with, if provided, they will override the attributes from the data

        See Also
        --------
        convert.to_compound 
        
        Spectrum

        """
        
        # define the attributes of the class
        self.id = None
        self.spectrum = None
        self.structure = None
        self.is_known = None
        self.name = None
        self.peak_fragments_map = None
        self.distances = None
        self.usi = None
        self.adduct_mass = None
        self.additional_attributes = {}
        self.exact_mass = None

        if incoming_data is None and len(kwargs) == 0:
            return

        # attempt to initialize the class with the provided data
        if incoming_data is not None:
            convert.to_compound(incoming_data, use_object = self, **kwargs)
        
        # update the attributes with the provided keyword arguments
        self.update(**kwargs)

        # TODO: write setters for spectrum, structure to warn the user to update the dependent attributes

    
    def clear(self):
        """Clear the compound data and reset all the attributes to None"""
        self.id = None
        self.spectrum = None
        self.structure = None
        self.is_known = None
        self.name = None
        self.peak_fragments_map = None
        self.distances = None
        self.usi = None
        self.adduct_mass = None
        self.additional_attributes = {}
        
        # # remove any additional attributes
        # all_keys = list(self.__dict__.keys())
        # for key in all_keys:
        #     if key not in ['id', 'spectrum', 'structure', 'is_known', 'name', 'peak_fragments_map', 'distances', 'usi', 'adduct_mass']:
        #         delattr(self, key)
        
    

    def update(self, structure = None, id: str = None, spectrum: Spectrum = None, usi: str = None, 
               adduct_mass: float = None, is_known: bool = None, name: str = None,
               peak_fragments_map: dict = None, distances: dict = None, **kwargs):
        """Update the attributes of the class

        Parameters
        ----------
        structure (Chem.Mol, Smiles, InChi): The structure of the compound
        id (str): The id of the compound
        spectrum (Spectrum): an instance of Spectrum containing peak information [(mz, intensity)], precursor mass, charge, and adduct for mass spectrumetry data
        usi (str): The USI of the compound
        adduct_mass (float): The adduct mass, derived from the Adduct
        is_known (bool): A boolean indicating whether the compound is known, if set to False, annotators or other parts of the code will treat this compound as unknown, if not provided but the structure is provided, it will be set to True
        name (str): The name of the compound
        peak_fragments_map (dict): A dictionary mapping peaks to fragments
        distances (dict): A dictionary of distances between every pair of atoms in the compound, if not provided, it will be calculated from the structure
        **kwargs: Additional data
            - A use case is for the scenarios where the *data* parameter is provided, these arguments will be used to parse and clean the data
        """
        # convert keys to lowercase
        lower_kwargs = {key.lower(): value for key, value in kwargs.items()}
        try:
            temp_structure = _get_molecule(structure, **lower_kwargs)
            self.structure = temp_structure if temp_structure is not None else self.structure
        except Exception:
            pass
        self.id = id if id is not None else self.id
        if spectrum is not None or "precursor_mz" in lower_kwargs:
            try:
                spectrum = convert.to_spectrum(spectrum, **lower_kwargs)
            except Exception:
                spectrum = None
        
        if spectrum is not None and spectrum.mz is not None:
            self.spectrum = spectrum
        self.usi = usi if usi is not None else self.usi
        self.adduct_mass = adduct_mass if adduct_mass is not None else self.adduct_mass
        self.is_known = is_known if is_known is not None else self.is_known
        self.name = name if name is not None else self.name
        if self.name is None:
            self.name = lower_kwargs.get('compound_name', None)
        self.peak_fragments_map = peak_fragments_map if peak_fragments_map is not None else self.peak_fragments_map
        self.distances = distances if distances is not None else self.distances

        # update the rest of the attributes
        self.additional_attributes.update(kwargs)
        
        self._parse_data()
    

    def _parse_data(self):
        """ Parse missing and verify the data of the class"""
        if self.is_known is None:
            self.is_known = (self.structure is not None)
            
        if self.structure is not None:
            self.exact_mass = rdMolDescriptors.CalcExactMolWt(self.structure)
        
        if self.distances is None and self.structure is not None:
            Chem.rdmolops.GetDistanceMatrix(self.structure)
        
        if self.adduct_mass is None and self.spectrum is not None:
            self.adduct_mass = general_utils.get_adduct_mass(self.spectrum.adduct)
        
        if self.is_known is None and self.structure is not None:
            self.is_known = True
        
        # TODO: check for a valid compound
        # what is needed for a compound:
        # id, spectrum.peaks, spectrum.precursor_mass, spectrum.charge, spectrum.adduct
    

    def __str__(self):
        valuable_data_keys = ['id', 'name', 'usi']
        strings = [f"{key}: {self.__dict__[key]}" for key in valuable_data_keys if self.__dict__[key] is not None]
        joined = ', '.join(strings)

        result = f"Compound({joined}) with {len(self.spectrum.mz)} peaks"
        if self.structure is not None:
            result += f" and structure {Chem.MolToSmiles(self.structure)}"

        return result
    

    def copy(self):
        """Return a copy of the compound"""
        copied_compound = Compound()
        convert.to_compound(self, use_object = copied_compound)
        return copied_compound
    

    def get_meta_data(self):
        """ Get the meta data of the compound

        Returns:
            dict: A dictionary containing the meta data of the compound
        """
        description = {
            "num_peaks": len(self.spectrum.mz),
            "adduct": self.spectrum.adduct,
            "precursor_mz": self.spectrum.precursor_mass,
            "charge": self.spectrum.charge,
        }

        if self.name is not None:
            description["name"] = self.name
        
        if self.is_known:
            description['num_aromatic_rings'] = rdMolDescriptors.CalcNumAromaticRings(self.structure)
            description['num_atoms'] = self.structure.GetNumAtoms()
            description['num_bonds'] = self.structure.GetNumBonds()
            description['num_rings'] = rdMolDescriptors.CalcNumRings(self.structure)

            # TODO: add peak annotation information
        
        return description
    
    
    def print_peaks_to_fragments_map(self, peaks: list = None):
        """Print the peaks to fragments map
        
        Parameters
        ----------
        peaks : list, optional (default: None)
            The list of peaks to print the fragments for, if None, print all peaks
        """
        if peaks is None:
            peaks = range(len(self.spectrum.mz))
        
        for peak in peaks:
            print(f"Peak {peak}: {self.spectrum.mz[peak]}, Fragments: {self.peak_fragments_map[peak]}")
        print()
        
    
    def find_existance(self, peakids: list):
        """
        For each atom, and for each peak in the list, find the fragments that the atom is part of
        
        Parameters
        ----------
            peakids : list of peak ids
        
        Returns
        -------
            existance : list of dicts for each atom, each dict contains the peak ids as keys and the fragments as values
            
        """
        
        existance = [dict() for i in range(len(self.structure.GetAtoms()))]
        for peak in peakids:
            for fragment in self.peak_fragments_map[peak]:
                # get all the bits that are 1 in the fragment
                bin_fragment = bin(fragment)
                len_fragment = len(bin_fragment)
                hitAtoms = [len_fragment-i-1 for i in range(len(bin_fragment)) if bin_fragment[i] == '1']
                for atom in hitAtoms:
                    if peak not in existance[atom]:
                        existance[atom][peak] = []
                    existance[atom][peak].append(fragment)
        return existance
    
    
    def calculate_contribution_atom_in_peak(self, atom: int, peakindx:int, existance_data:list, CI:bool = False, CPA:bool = True, CFA:bool = True):
        """Calculates the contribution of an atom to a peak

        Parameters
        ----------
        atom : int
            The index of the atom
        peakindx : int
            The index of the peak
        existance_data : list of dicts
            The existance data for the atoms
        CI : bool, optional (default: False)
            Calculate the intensity factor
        CPA : bool, optional (default: True)
            Calculate the atom peak ambiguity factor
        CFA : bool, optional (default: True)
            Calculate the fragment ambiguity factor
        """
        contribution = 0
        if peakindx not in existance_data[atom]:
            return contribution
        
        intensity_factor = 1
        atom_peak_ambiguity_factor = 1
        fragment_ambiguity_factor = 1

        if CI:
            intensity_factor = self.spectrum.intensity[peakindx]
        if CPA:
            atom_peak_ambiguity_factor = 1/len(self.peak_fragments_map[peakindx])

        for frag in existance_data[atom][peakindx]:
            if CFA:
                fragment_ambiguity_factor = 1/(bin(frag).count("1"))
            
            contribution += intensity_factor * atom_peak_ambiguity_factor * fragment_ambiguity_factor
        
        return contribution
    
    
    def calculate_contributions(self, peakids, CI = False, CPA = True, CFA = True, CPE = True):
        """ 
        input:
        peakids: list of peak ids
        CI: (Consider_Intensity) bool, if True, the intensity of the peaks is considered (default: False)
        CPA: (Consider_Peak_Ambiguity) bool, if True, the peak ambiguity (number of fragments assigned to a peak) is considered (default: True)
        CFA: (Consider_Fragment_Ambiguity) bool, if True, the fragment ambiguity (number of atoms in fragment) is considered (default: True)
        CPA: (Consider_Peak_Entropy) bool, if True, the peak entropy (how ambiguis the fragments are) is considered (default: True
        """
        num_atoms = len(self.structure.GetAtoms())
        existance_data = self.find_existance(peakids)
        contributions = [0 for i in range(num_atoms)]
        peak_atom_contributions = np.zeros((len(peakids), num_atoms))
        for i, peak in enumerate(peakids):
            for atom in range(num_atoms):
                peak_atom_contributions[i][atom] = self.calculate_contribution_atom_in_peak(atom, peak, existance_data, CI=CI, CPA=CPA, CFA=CFA)
        
        if CPE:
            peak_entropies = np.zeros(len(peakids))
            for i in range(len(peakids)):
                peak_entropies[i] = 1 - general_utils.entropy(peak_atom_contributions[i])
            
            # peak_entropies = peak_entropies / np.max(peak_entropies)
        else:
            peak_entropies = np.ones(len(peakids))
            
        
        for i in range(num_atoms):
            for j in range(len(peakids)):
                contributions[i] += peak_atom_contributions[j][i] * peak_entropies[j]
        
        return contributions


    def filter_fragments_by_atoms(self, atoms: list, peaks: list = None):
        """Filter the fragments by the atoms, remove fragments that do not contain at least one of the atoms
        Parameters
        ----------
            atoms: a list of atoms to filter the fragments
            peaks: a list of peaks to filter their fragments, if None, use all peaks
        
        Returns
        -------
            updated: the number of updated peaks
        """
        if peaks is None:
            peaks = [i for i in range(len(self.peaks))]
        
        updated = 0
        for i in peaks:
            updated_fragments = set()
            for fragment in self.peak_fragments_map[i]:
                for atom in atoms:
                    if 1 << atom & fragment:
                        updated_fragments.add(fragment)
                        break

            if len(updated_fragments) != len(self.peak_fragments_map[i]):
                updated += 1
            
            self.peak_fragments_map[i] = updated_fragments

        return updated
    
    
    def calculate_peak_annotation_ambiguity(self, peaks: list=None):
        """Calculate the peak annotation ambiguity
        
        Parameters
        ----------
            peaks : list
                a list of peaks to calculate the ambiguity for, if None, use all peaks
        Returns
        -------
            (float, float) : a tuple of two values (ambiguity, ratio)
                ambiguity: the average number of fragments per annotated peaks
                ratio: the ratio of annotated peaks to all peaks
        """
        if peaks is None:
            peaks = [i for i in range(len(self.peaks))]
        ambiguity = 0
        annotated_peaks = 0
        for peak in peaks:
            if len(self.peak_fragments_map[peak]) > 0:
                annotated_peaks += 1
                ambiguity += len(self.peak_fragments_map[peak])
        if annotated_peaks == 0:
            return -1, 0
        if len(peaks) == 0:
            return -1, -1
        return ambiguity / annotated_peaks, annotated_peaks / len(peaks)


    def calculate_annotation_entropy(self, peaks: list=None):
        """Calculate the entropy of the annotation
        
        Parameters
        ----------
            peaks : list
                a list of peaks to calculate the entropy for, if None, use all peaks
        
        Returns
        -------
            float : the entropy of the annotation
        """
        if peaks is None:
            peaks = [i for i in range(len(self.peaks))]
        peak_entropies = [0 for i in range(len(self.peaks))]
        n = len(self.structure.GetAtoms())
        for peak in peaks:
            atoms_appearance = [0 for i in range(n)]
            for fragment in self.peak_fragments_map[peak]:
                for atom in range(n):
                    if 1 << atom & fragment:
                        atoms_appearance[atom] += 1
            entropy = 0
            for atom in range(n):
                if atoms_appearance[atom] > 0:
                    p = atoms_appearance[atom] / len(self.peak_fragments_map[peak])
                    entropy -= p * math.log(p)
            peak_entropies[peak] = entropy
        if len(peak_entropies) == 0:
            return -1
        entropy = sum(peak_entropies) / len(peak_entropies)
        return entropy
