from modifinder.engines.Abtracts import AnnotationEngine
from modifinder.classes.Compound import Compound
from modifinder.classes.EdgeDetail import EdgeDetail, MatchType
from typing import List, Tuple
from modifinder.engines.annotation.magma import (
    fragmentation_py,
    rdkit_engine as rdkit_engine,
)
from modifinder.utilities.general_utils import mims, is_submolecule
from modifinder.utilities.mol_utils import get_modification_nodes
import networkx as nx
from rdkit import Chem
from msbuddy import assign_subformula


class MAGMaAnnotationEngine(AnnotationEngine):
    def __init__(self, breaks = 2, max_water_losses = 2, ionisation_mode = 0, mz_tolerance = 0.1, ppm_tolerance = 40, **kwargs):
        """
        Initializes the MAGMa annotation engine

        Parameters:
            :kwargs: additional arguments
        """
        self.breaks = breaks
        self.max_water_losses = max_water_losses
        self.ionisation_mode = ionisation_mode
        self.mz_tolerance = mz_tolerance
        self.ppm_tolerance = ppm_tolerance
        self.args = kwargs


    def annotate(self, network: nx.DiGraph, annotate_all: bool = False, **kwargs):
        """
        Annotates the network using the MAGMa model

        Parameters:
            :network (Network): the network to be annotated
            :annotate_all (bool): whether to annotate all compounds or only the ones that have not been annotated yet
            :kwargs: additional arguments
        """
        for node in network.nodes:
            compound = network.nodes[node]["compound"]
            if compound is not None and compound.is_known:
                if annotate_all or (compound.spectrum.peak_fragment_dict is None or len(compound.spectrum.peak_fragment_dict) == 0):
                    self.annotate_single(compound, modify_compound=True, **kwargs)
        
        # refine by helpers
        for edge in network.edges:
            edge_detail = network.edges[edge]["edgeDetail"]
            if edge_detail is not None:
                first_compound = network.nodes[edge[0]]["compound"]
                second_compound = network.nodes[edge[1]]["compound"]
                if first_compound.is_known and second_compound.is_known:
                    try:
                        self.refine_annotations_by_helper(network.nodes[edge[0]]["compound"], network.nodes[edge[1]]["compound"], edge_detail, modify_compound = True)
                        self.refine_annotations_by_helper(network.nodes[edge[1]]["compound"], network.nodes[edge[0]]["compound"], edge_detail, modify_compound = True)
                    
                    except Exception as e:
                        raise e


    def annotate_single(
        self, compound: Compound, modify_compound: bool = True, **kwargs
    ) -> List[Tuple[int, List[str]]]:
        """
        provides annotation for the peaks in a single compound

        Parameters:
            :compound (Compound): the compound to be annotated
            :modify_compound (bool): whether to modify the passed compound with the annotations
            :kwargs: additional arguments

        Returns:
            :Mapping[int, List[str]]: a dictionary with the indices of the peaks as keys and the list of annotations as values
        """

        if "breaks" not in kwargs:
            kwargs["breaks"] = self.breaks
        if "mz_tolerance" not in kwargs:
            kwargs["mz_tolerance"] = self.mz_tolerance
        if "ppm_tolerance" not in kwargs:
            kwargs["ppm_tolerance"] = self.ppm_tolerance

        for item in ["breaks", "mz_tolerance", "ppm_tolerance"]:
            if item not in kwargs or kwargs[item] is None:
                raise ValueError(
                    f"Missing parameters for MAGMa annotation engine {item}"
                )

        fragmentation_instance = fragmentation_py.FragmentEngine(
            compound.structure, kwargs["breaks"], max_water_losses=self.max_water_losses, ionisation_mode=self.ionisation_mode, skip_fragmentation=False, molcharge=compound.spectrum.precursor_charge
        )
        fragmentation_instance.generate_fragments()

        # generate peak to fragment map
        base_precision = 1 + kwargs["ppm_tolerance"] / 1000000
        peak_fragment_dict = {int(key): set() for key in compound.spectrum.mz_key}
        for mz_key in compound.spectrum.mz_key:
            search_weight = (mz_key / 1e6) - compound.spectrum.adduct_mass
            annotations = fragmentation_instance.find_fragments(
                search_weight, 0.1, base_precision, kwargs["mz_tolerance"]
            )
            for annotation in annotations:
                try:
                    peak_fragment_dict[int(mz_key)].add(annotation[0])
                except KeyError:
                    # Get keys within 0.002 m/z for debugging
                    print(f"KeyError for m/z {mz_key}. Available keys: {[key for key in peak_fragment_dict.keys() if abs(key - mz_key) <= 1e4]}", flush=True)
                    print(f"Sample of some keys: {list(peak_fragment_dict.keys())[:10]}", flush=True)
                    raise KeyError(f"m/z {mz_key} not found in peak_fragment_dict")

        # refine by formula
        # peak_fragment_dict = self.refine_annotations_by_formula(compound, peak_fragment_dict, modify_compound = False)
        
        if modify_compound:
            compound.spectrum.peak_fragment_dict = peak_fragment_dict

        return peak_fragment_dict


    def get_fragment_info(self, Compound: Compound, fragment):
        """
        converts a fragment to a SMILES string

        Parameters
        ----------
            Compound (Compound): the compound
            fragment: the fragment representation

        Returns
        -------
            tuple(atomlist -> List[int], edge_list -> List[Tuple[int, int]], formula -> str, smiles -> str): a tuple containing the atom list, the edge list, the formula and the SMILES string
        """
        
        atomstring = ""
        atomlist = []
        edgeList = []
        elements = dict([(e, 0) for e in mims.keys()])
        natoms = Compound.structure.GetNumAtoms()
        for atom in range(natoms):
            if 1 << atom & fragment:
                atomstring += "," + str(atom)
                atomlist.append(atom)
                elements[rdkit_engine.GetAtomSymbol(Compound.structure, atom)] += 1
                elements["H"] += rdkit_engine.GetAtomHs(Compound.structure, atom)

        for bond in Compound.structure.GetBonds():
            if bond.GetBeginAtomIdx() in atomlist and bond.GetEndAtomIdx() in atomlist:
                edgeList.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))

        formula = ""
        for el in mims.keys():
            nel = elements[el]
            if nel > 0:
                formula += el
            if nel > 1:
                formula += str(nel)

        return (
            atomlist,
            edgeList,
            formula,
            rdkit_engine.FragmentToSmiles(Compound.structure, atomlist),
        )
    
    
    def refine_annotations_by_formula(self, compound: Compound, peak_fragment_dict, modify_compound: bool = True):
        """
        Refines the annotations of a compound by using the formula from msbuddy

        See Also
        --------
        `MSBuddy <https://github.com/Philipbear/msbuddy>`_
        
        
        Parameters
        ----------
            compound (Compound): the compound
            
        """
        raise NotImplementedError("This function has not been tested in the new verion")
        print("Refining the annotations by formula", flush=True)
        
        structure = compound.structure
        spectrum = compound.spectrum
        if structure is None:
            return peak_fragment_dict
        
        main_compound_formula = Chem.rdMolDescriptors.CalcMolFormula(structure)
        peak_mz = spectrum.mz_key / 1e6
        
        if len(spectrum.mz_key) == 0:
            return peak_fragment_dict
        
        subformla_list = assign_subformula(peak_mz,
                                        precursor_formula=main_compound_formula, adduct=spectrum.adduct,
                                        ms2_tol=self.ppm_tolerance, ppm=True, dbe_cutoff=-1.0)
        
        new_peak_fragment_dict = {int(key): set() for key in spectrum.mz_key_ids}
        for idx, key in enumerate(spectrum.mz_key):
            possibilites = set()
            for subformula in subformla_list[idx].subform_list:
                formula = subformula.formula
                # formula = utils.remove_adduct_from_formula(formula, self.Adduct)
                
                # find the fragments that contains the formula
                for frag_id in peak_fragment_dict[int(key)]:
                    molSubFormula = self.get_fragment_info(compound, frag_id)[2]
                    if is_submolecule(molSubFormula, formula, self.args["formula_ignore_H"]) and is_submolecule(formula, molSubFormula, self.args["formula_ignore_H"]):
                        possibilites.add(frag_id)
            
            if len(possibilites) > 0:
                new_peak_fragment_dict[int(key)] = possibilites
        

        print(f"Refining the annotations by formula. Starting with {sum([len(peak_fragment_dict[i]) for i in peak_fragment_dict.keys()])} annotations and after refining, {sum([len(new_peak_fragment_dict[i]) for i in new_peak_fragment_dict.keys()])} annotations remain", flush=True)

        if modify_compound:
            compound.spectrum.peak_fragment_dict = new_peak_fragment_dict
        
        return new_peak_fragment_dict
    
    
    def refine_annotations_by_helper(self, main_compound: Compound, helper_compound: Compound, edgeDetail: EdgeDetail, modify_compound: bool = True):
        """
        Refines the annotations of a compound by using the helper compound
        """
        if not modify_compound:
            main_compound = main_compound.copy()
            
        if helper_compound.exact_mass < main_compound.exact_mass:
            modification_site = get_modification_nodes(main_compound.structure, helper_compound.structure, True)
            shifted_peaks = [match.second_peak_mz for match in edgeDetail.matches if match.match_type == MatchType.shifted]
            sub_match_indices = main_compound.structure.GetSubstructMatch(helper_compound.structure)
            mapping = dict()
            for i, atom in enumerate(sub_match_indices):
                mapping[i] = atom
            unshifted = [(match.second_peak_mz, match.first_peak_mz) for match in edgeDetail.matches if match.match_type == MatchType.unshifted]
        else:
            modification_site = get_modification_nodes(helper_compound.structure, main_compound.structure, False)
            shifted_peaks = [match.first_peak_mz for match in edgeDetail.matches if match.match_type == MatchType.shifted]
            sub_match_indices = helper_compound.structure.GetSubstructMatch(main_compound.structure)
            mapping = dict()
            for i, atom in enumerate(sub_match_indices):
                mapping[atom] = i
            unshifted = [(match.first_peak_mz, match.second_peak_mz) for match in edgeDetail.matches if match.match_type == MatchType.unshifted]
        
        main_compound.filter_fragments_by_atoms_and_peak_ids(modification_site, shifted_peaks)
        
        # filter unshifted peaks by intersection
        
        for edge in unshifted:
            helper_fragments = set()
            for fragment in helper_compound.spectrum.peak_fragment_dict[int(edge[1])]:
                new_fragment = 0
                for i in range(len(helper_compound.structure.GetAtoms())):
                    if 1 << i & fragment:
                        if i not in mapping:
                            new_fragment = -1
                            break
                        else:
                            new_fragment += 1 << mapping[i]
                if new_fragment != -1:
                    helper_fragments.add(new_fragment)
            
            fragments_assigned_to_peak = set(main_compound.spectrum.peak_fragment_dict[int(edge[0])])
            main_compound.spectrum.peak_fragment_dict[int(edge[0])] = main_compound.spectrum.peak_fragment_dict[int(edge[0])].intersection(helper_fragments)
             
        return main_compound.spectrum.peak_fragment_dict
        
        
            
        
        
        
        
        
        

        