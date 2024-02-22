from rdkit import Chem
from msbuddy import assign_subformula
from . import fragmentation_py as fragmentation_py
import copy
import json
from . import utils_n as utils
from . import handle_network as handle_network
import re
import math

important_arguments = ["peaks", "Adduct", "Precursor_MZ", "Charge"]

default_args = {
    "ppm": (int, 40), # the ppm tolerance used for matching peaks
    "mz_tolerance": (float, 1), # the mz tolerance for matching peaks
    "filter_peaks_method": (str, "both"), # can be "intensity", "top_k", "none", "both"
    "filter_peaks_variable": (float, 0.01), # if intensity, the percentage of the highest peak, if top_k, the k (int), if both, the percentage of the highest peak
    "fragmentation_depth": (int, 2), # -1 means auto, otherwise, the number of breaks
    "formula_ignore_H": (bool, True), # whether to ignore H when comparing formulas
    "should_fragment": (bool, True), # whether to fragment the compound
}


class Compound:
    """A class to represent a compound."""

    def __init__(self, data, structure=None, args={}):

        """Initialize the compound."""
        
        # Adding the args -----------------------------------------------
        self.args = {}
        ## set the default args and update it with the args passed in
        for arg in default_args:
            self.args[arg] = default_args[arg][1]
        self.args.update(args)
        for arg in list(args.keys()):
            # if its a string and can be converted to float, convert it
            if type(args[arg]) == str:
                try:
                    self.args[arg] = float(args[arg])
                except:
                    self.args[arg] = args[arg]

        # defining the attributes ----------------------------------------
        self.name = None
        self.accession = None
        self.library_membership = None
        self.peaks = None
        self.Adduct = None
        self.Precursor_MZ = None
        self.Charge = None
        self.fragments = None
        self.numFrag = None
        self.peak_fragments_map = None
        self.Smiles = None
        self.structure = None
        self.distances = None
        
        # if data is a string (USI or accession), get the data from the network
        if type(data) == str:
            try:
                # if data contains substring accession
                if data.count(':') == 0 or "accession" in data:
                    accession = data
                    if "accession" in data:
                        accession = data.split(':')[-1]
                    self.accession = accession
                    data = handle_network.getDataFromAccession(accession)
                    data_part_one = data["annotations"][0]
                    data_part_two = data["spectruminfo"]
                    data = {**data_part_one, **data_part_two}
                else:
                    data = handle_network.getDataFromUsi(data)
            except:
                raise ValueError("Invalid input, error getting data from network")
        
        # setting the attributes -----------------------------------------
        if "Compound_Name" in data:
            self.name = data["Compound_Name"]
        if "library_membership" in data:
            self.library_membership = data["library_membership"]
        if "spectrum_id" in data:
            self.accession = data["spectrum_id"]
            
        for arg in important_arguments:
            if not arg in data and not arg.lower() in data:
                if arg != "peaks":
                    raise ValueError("Missing argument: " + arg)
                elif "peaks_json" in data:
                    self.peaks = json.loads(data["peaks_json"])
                else:
                    raise ValueError("Missing argument: " + arg)
            else:
                if arg.lower() in data:
                    setattr(self, arg, data[arg.lower()])
                else:
                    setattr(self, arg, data[arg])

        # adjusting the attributes ---------------------------------------
        self.Charge = int(self.Charge)
        self.Precursor_MZ = float(self.Precursor_MZ)
        self.Adduct = utils.parse_adduct(self.Adduct)

        self.peaks = utils.filter_peaks(
            self.peaks,
            self.args["filter_peaks_method"],
            self.args["filter_peaks_variable"],
            self.Precursor_MZ,
            self.Charge,
        )


        # set the smiles and structure-----------------------------------
        if "Smiles" in data or (structure != None and type(structure) == str and len(structure) > 0):
            if "Smiles" in data:
                self.Smiles = data["Smiles"]
            else:
                self.Smiles = structure

        self.structure = None
        if structure != None and type(structure) != str:
            self.structure = structure
        elif self.Smiles != None:
            self.structure = Chem.MolFromSmiles(self.Smiles)

        # perform fragmentation------------------------------------------
        if self.structure != None:
            self.distances = Chem.rdmolops.GetDistanceMatrix(self.structure)
            
            if self.args["should_fragment"]:
                self.generate_fragments()
    
    def generate_fragments(self):
        if self.args["fragmentation_depth"] == -1:
            breaks = 4
            if (self.structure.GetNumAtoms() < 20):
                breaks = 2
            elif (self.structure.GetNumAtoms() < 40):
                breaks = 3
        else:
            breaks = self.args["fragmentation_depth"]
            # if string, convert to int
            if type(breaks) != int:
                breaks = int(breaks)
            if (self.structure.GetNumAtoms() > 80):
                breaks = min(breaks, 2)
        self.fragments = fragmentation_py.FragmentEngine(
            Chem.MolToMolBlock(self.structure), breaks, 2, 1, 0, 0
        )
        self.numFrag = self.fragments.generate_fragments()
        self.generate_peak_to_fragment_map()


    def generate_peak_to_fragment_map(self):
        base_precision = 1 + self.args["ppm"] / 1000000
        self.peak_fragments_map = [set() for i in range(len(self.peaks))]
        for i in range(len(self.peaks)):
            annotations = self.fragments.find_fragments(
                self.peaks[i][0], 0.1, base_precision, self.args["mz_tolerance"]
            )
            for annotation in annotations:
                self.peak_fragments_map[i].add(annotation[0])

    def apply_sirius(self, sirius, add_adduct=True):
        """Apply the sirius results to the compound and filter the possible annotations"""
        for i, peak in enumerate(self.peaks):
            index = utils.find_mz_in_sirius(
                sirius["fragments"], peak[0], self.args["mz_tolerance"], self.args["ppm"]
            )
            if index == -1:
                continue
            helper_peak_formula = sirius["fragments"][index]["molecularFormula"]
            if add_adduct:
                helper_peak_formula = utils.add_adduct_to_formula(helper_peak_formula, self.Adduct)
            possibilites = set()
            for frag_id in self.peak_fragments_map[i]:
                molSubFormula = self.fragments.get_fragment_info(frag_id, 0)[2]
                if utils.is_submolecule(molSubFormula, helper_peak_formula, self.args["formula_ignore_H"]) and utils.is_submolecule(helper_peak_formula, molSubFormula, self.args["formula_ignore_H"]):
                    possibilites.add(frag_id)
            
            if len(possibilites) > 0:
                self.peak_fragments_map[i] = possibilites

    def apply_filter(self, filter_method, filter_value):
        self.peaks = utils.filter_peaks(
            self.peaks, filter_method, filter_value, self.Precursor_MZ, self.Charge
        )

    def remove_large_peaks(self, weight=None):
        if weight == None:
            weight = self.Precursor_MZ - 1.6

        peaks = []
        for i in range(len(self.peaks)):
            if self.peaks[i][0] < weight:
                peaks.append(self.peaks[i])
        self.peaks = peaks

    def get_peak_index(self, peak_weight, args={}):
        args.update(self.args)
        ind = []
        for i in range(len(self.peaks)):
            diff = abs(self.peaks[i][0] - peak_weight)
            if diff < args["mz_tolerance"] and ((diff / peak_weight) * 1000000) < args["ppm"]:
                ind.append(i)
        return ind


    def filter_fragments_by_atoms(self, atoms, peaks):
        """Filter the fragments by the atoms, remove fragments that do not contain at least one of the atoms
        Args:
            atoms: a list of atoms to filter the fragments
            peaks: a list of peaks to filter their fragments, if None, use all peaks
        """
        if peaks == None:
            peaks = [i for i in range(len(self.peaks))]
        for i in peaks:
            updated_fragments = set()
            for fragment in self.peak_fragments_map[i]:
                for atom in atoms:
                    if 1 << atom & fragment:
                        updated_fragments.add(fragment)
                        break
            self.peak_fragments_map[i] = updated_fragments


    def vote_for_fragments(self, peak_index, fragments, voter="none", extra_args={}):
        """Vote for the fragments, update the peak_fragments_map with the intersection of the fragments, and extend the intersection to other peaks if extend_to_other_peaks is true"""
        
        print("here in vote for fragments", peak_index, fragments)
        if peak_index == -1:
            print("in vote for fragment, Peak not found")
            return

        # get the intersection of the fragments
        intersection = self.peak_fragments_map[peak_index].intersection(fragments)
        self.peak_fragments_map[peak_index] = intersection
        extra_args.update(self.args)
        # if kwargs has extend_to_other_peaks, extend the intersection to other peaks
        if "extend_to_other_peaks" in extra_args and extra_args["extend_to_other_peaks"]:
            # find atoms in the intersection
            intersect_atoms = set()
            for atom in range(len(self.structure.GetAtoms())):
                exists = True
                for fragment in intersection:
                    if not 1 << atom & fragment:
                        exists = False
                        break
                if exists:
                    intersect_atoms.add(atom)
            
            self.filter_fragments_by_atoms(intersect_atoms, None)
        return
    

    def calculate_peak_annotation_ambiguity(self, peaks=None):
        """Calculate the peak annotation ambiguity
        Args:
            peaks: a list of peaks to calculate the ambiguity for, if None, use all peaks
        Returns:
            ambiguity: the average number of fragments per annotated peaks
            ratio: the ratio of annotated peaks to all peaks
        """
        if peaks == None:
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
    
    def apply_msbuddy(self):
        if not "M+H" in self.Adduct:
            raise ValueError("Adduct not supported")
        main_compound_formula = Chem.rdMolDescriptors.CalcMolFormula(self.structure)
        peak_mz = [peak[0] for peak in self.peaks]
        if len(peak_mz) == 0:
            return
        subformla_list = assign_subformula(peak_mz,
                                        precursor_formula=main_compound_formula, adduct="[M+H]+",
                                        ms2_tol=self.args["ppm"], ppm=True, dbe_cutoff=-1.0)
        
        for i in range(len(self.peaks)):
            possibilites = set()
            for subformula in subformla_list[i].subform_list:
                formula = subformula.formula
                # formula = utils.remove_adduct_from_formula(formula, self.Adduct)
                
                # find the fragments that contains the formula
                for frag_id in self.peak_fragments_map[i]:
                    molSubFormula = self.fragments.get_fragment_info(frag_id, 0)[2]
                    if utils.is_submolecule(molSubFormula, formula, self.args["formula_ignore_H"]) and utils.is_submolecule(formula, molSubFormula, self.args["formula_ignore_H"]):
                        possibilites.add(frag_id)
            
            if len(possibilites) > 0:
                self.peak_fragments_map[i] = possibilites
    
    def calculate_annotation_entropy(self, peaks=None):
        """Calculate the entropy of the annotation
        Args:
            peaks: a list of peaks to calculate the entropy for, if None, use all peaks
        Returns:
            entropy: the entropy of the annotation
        """
        if peaks == None:
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

    def apply_helper(self, helper_compound, shifted, unshifted, unshifted_mode = "union"):
        """use helper compound to update the peak_fragments_map
        Args:
            helper_compound: the helper compound
            shifted: array of pairs of the shifted peaks, each pair is (self_peak_index, helper_peak_index)
            unshifted: array pf pairs of the unshifted peaks, each pair is (self_peak_index, helper_peak_index)
            unshifted_mode: the mode to update the unshifted peaks, can be "union", "intersection", None
        """
        if not (self.structure.HasSubstructMatch(helper_compound.structure) or helper_compound.structure.HasSubstructMatch(self.structure)):
            print("helper compound and main compound do not have common substructure, no change applied")
            return


        if helper_compound.Precursor_MZ < self.Precursor_MZ:
            modification_site = utils.calculateModificationSites(self.structure, helper_compound.structure, True)
        else:
            modification_site = utils.calculateModificationSites(helper_compound.structure, self.structure, False)
        
        # update the shifted peaks
        ## get shifted indices of self
        shifted_indices = [_[0] for _ in shifted]
        self.filter_fragments_by_atoms(modification_site, shifted_indices)

        if unshifted_mode == None or unshifted_mode == "none":
            return
        # update the unshifted peaks
        ## get a mapping from the helper atoms index to the self aoms index
        if helper_compound.Precursor_MZ < self.Precursor_MZ:
            sub_match_indices = self.structure.GetSubstructMatch(helper_compound.structure)
            mapping = {}
            for i, atom in enumerate(sub_match_indices):
                mapping[i] = atom
        else:
            sub_match_indices = helper_compound.structure.GetSubstructMatch(self.structure)
            mapping = {}
            for i, atom in enumerate(sub_match_indices):
                mapping[atom] = i

        for peak in unshifted:
            helper_peak_fragment_map = set()
            for fragment in helper_compound.peak_fragments_map[peak[1]]:
                new_fragment = 0
                for i in range(len(helper_compound.structure.GetAtoms())):
                    if 1 << i & fragment:
                        if i not in mapping:
                            new_fragment = -1
                            break
                        else:
                            new_fragment += 1 << mapping[i]
                if new_fragment != -1:
                    helper_peak_fragment_map.add(new_fragment)
            if unshifted_mode == "union":
                self.peak_fragments_map[peak[0]] = self.peak_fragments_map[peak[0]].union(helper_peak_fragment_map)
            elif unshifted_mode == "intersection":
                self.peak_fragments_map[peak[0]] = self.peak_fragments_map[peak[0]].intersection(helper_peak_fragment_map)
            else:
                raise ValueError("unshifted_mode not supported")
        return


    def apply_iceberg(self, iceberg):
        # TODO: implement
        pass
