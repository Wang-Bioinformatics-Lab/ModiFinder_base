"""Base class for ModiFinder.

This class is used to create a ModiFinder object.
The object can be used to get information about unknown compounds in the network using the known compounds.

This class Builds and maintains a network of compounds where the nodes are compounds (known and unknown)
"""

from typing import List

from .. import convert as convert
from modifinder.classes.Compound import Compound
from modifinder.classes.EdgeDetail import EdgeDetail, MatchType
from modifinder.engines import AlignmentEngine, AnnotationEngine
from modifinder.engines.annotation.MAGMaAnnotationEngine import MAGMaAnnotationEngine
from modifinder.engines.alignment.CosineAlignmentEngine import CosineAlignmentEngine
from modifinder.utilities import visualizer as mf_vis

from modifinder.exceptions import (
    ModiFinderNotImplementedError,
    ModiFinderNotSolvableError,
)

import networkx as nx
import numpy as np
from modifinder.utilities import general_utils as general_utils


class ModiFinder:
    """
    Base class for compound network.

    ModiFinder Can be used in two scenarios
    1. between a known and an unknown compound. 
    2. in a network consisting of known and unknown compounds.

    In both scenarios, the ModiFinder object stores a directed graph where the nodes are compounds (known and unknown)
    and edges are directed (from smaller compound to the larger compound) holding the relationships between the pair (alignment, similarity, etc). The object also stores
    the unknown compounds in the network. The object can be used to get information about the unknown compound.

    Parameters
    ----------
    network : nx.DiGraph
        A networkx graph object with the nodes are identified by the compound ids and the "compound"
        attibute of a node is a Compound object. The edges are the relationships between the compounds
        and the "edgeDetail" attribute of an edge is an EdgeDetail object.

    unknowns : list
        A list of compound ids that are unknown in the network.

    See Also
    --------
    Compound
    EdgeDetail.EdgeDetail
    AlignmentEngine
    modifinder.AlignmentEngine
    modifinder.engines.AlignmentEngine
    modifinder.engines.Abstracts.AlignmentEngine

    Examples
    --------
    """

    def __init__(
        self,
        knownCompound: Compound = None,
        unknownCompound: Compound = None,
        edgeDetail: EdgeDetail = None,
        helpers: List[Compound] = None,
        network: nx.DiGraph = None,
        networkUnknowns: list = None,
        should_align: bool = True,
        alignmentEngine: AlignmentEngine = None,
        should_annotate: bool = True,
        annotationEngine: AnnotationEngine = None,
        ppm_tolerance: float = 40,
        **kwargs,
    ):
        """
        Initialize the ModiFinder object.

        **Use Case  1: a pair of known and unknown compound**

        **Use Case  2: a network of known and unknown compounds**

        Parameters
        ----------
        knownCompound : known compound (not optional in Use Case 1, ignored in Use Case 2)
            Data to create a known compound. The data can be a Compound object, or any data that can be converted to a Compound object.

        unknownCompound : unknown compound (not optional in Use Case 1, ignored in Use Case 2)
            Data to create an unknown compound. The data can be a Compound object, or any data that can be converted to a Compound object.

        edgeDetail : EdgeDetail object from smaller to larger compound (optional in Use Case 1, ignored in Use Case 2)
            the orientation of the match must be from the smaller compound to the larger compound.
            
        helpers : list of helper compounds for the known compound (optional in Use Case 1, ignored in Use Case 2)
            A list of helper compounds to help in the annotation refinement of the known compound.
        
        helpers_edgeDetails : list of EdgeDetail objects (optional in Use Case 1, ignored in Use Case 2)
            A list of edge details between the helper compounds and the known compound. The orientation of the match must be from the smaller compound to the larger compound.
            If not passed, the method will use the alignment engine to align the helper compounds with the known compound.
        
        network : nx.Graph (if passed, then Use Case 2 is used, if None, then Use Case 1 is used)
            A networkx graph object with the nodes are identified by the compound ids and the "compound"
            attibute of a node is a *Compound* object. The edges are the relationships between the compounds
            and the "edgeDetail" attribute of an edge is an *EdgeDetail* object.

        networkUnknowns : list (ignored in Use Case 1, optional in Use Case 2)
            A list of compound ids that are unknown in the network, if not passed, the unknowns will be driven from 'is_known' attribute
            of the compounds in the network.
        
        should_align : bool (optional, default True)
            If True, the method will align the network using the alignment engine.
        
        alignmentEngine : AlignmentEngine (optional)
            The alignment engine to use to align the unknown compound with the known compounds in the network.
            If not passed, the method will use the CosineAlignmentEngine.
        
        should_annotate : bool (optional, default True)
            If True, the method will annotate the network using the annotation engine.
            
        annotationEngine : AnnotationEngine (optional)
            The annotation engine to use to annotate the compounds in the network.
            If not passed, the method will use the MAGMaAnnotationEngine.
        
        kwargs : dict
            Additional parameters to pass to the alignment engine, the annotation engine, and other methods.


        See Also
        --------
        Compound
        EdgeDetail


        Examples
        --------

        """
        self.network = None
        self.knowns = None
        self.unknowns = None
        self.ppm_tolerance = ppm_tolerance

        # Add ppm_tolerance to kwargs
        kwargs["ppm_tolerance"] = ppm_tolerance
        
        if alignmentEngine is None:
            self.alignmentEngine = CosineAlignmentEngine(**kwargs)
        else:
            self.alignmentEngine = alignmentEngine
        
        if annotationEngine is None:
            self.annotationEngine = MAGMaAnnotationEngine(**kwargs)
        else:
            self.annotationEngine = annotationEngine

        if network:
            raise ModiFinderNotImplementedError("Use Case 2 is not fully implemented yet")
            self.network = network

            if networkUnknowns:
                self.unknowns = networkUnknowns
            else:
                self.unknowns = [
                    node
                    for node in network.nodes
                    if not network.nodes[node]["compound"].is_known
                ]

        else:
            self.network = nx.DiGraph()
            if unknownCompound is not None:
                unknownCompound = convert.to_compound(data=unknownCompound, **kwargs)
                unknownCompound.is_known = False
                self.network.add_node(unknownCompound.id, compound=unknownCompound)
                self.unknowns = [unknownCompound]
            
            if knownCompound is not None:
                knownCompound = convert.to_compound(data=knownCompound, **kwargs)
                self.network.add_node(knownCompound.id, compound=knownCompound)
                self.knowns = [knownCompound]
            
            
            if knownCompound is not None and unknownCompound is not None:
                self.add_adjusted_edge(knownCompound.id, unknownCompound.id, edgeDetail, **kwargs)

        
        if helpers is not None:
            for helper in helpers:
                if not isinstance(helper, Compound):
                    raise ValueError("Helper compounds must be of type Compound")
                try:
                    self.network.add_node(helper.id, compound=helper)
                    self.add_adjusted_edge(helper.id, knownCompound.id)
                except Exception as e:
                    print(f"Error adding helper compound: {e}")
                    raise e
        
        if should_align:
            self.re_align(self.alignmentEngine, **kwargs)
            
        if should_annotate:
            self.re_annotate(self.annotationEngine, **kwargs)

    def get_result(self):
        """
        Get the result of the ModiFinder object.

        The method will return a dictionary with the following keys:
        - probabilities: the probabilities of the atoms in the unknown compound.
        
        Returns
        -------
        dict
            A dictionary with the following keys:
            - probabilities: the probabilities of the atoms in the unknown compound.
        
        """

        unknownCompound = self.unknowns[0]
        knownCompound = self.knowns[0]

        if unknownCompound.structure is not None:
            if not (unknownCompound.structure.HasSubstructMatch(knownCompound.structure) or knownCompound.structure.HasSubstructMatch(unknownCompound.structure)):
                return None, None, None, None, "None of the structures are substructures of the other"
            if unknownCompound.structure.HasSubstructMatch(knownCompound.structure) and knownCompound.structure.HasSubstructMatch(unknownCompound.structure):
                return None, None, None, None, "Structures are the same"

        main_compound_peaks = [(knownCompound.spectrum.mz_key[i], knownCompound.spectrum.intensity[i]) for i in range(len(knownCompound.spectrum.mz_key))]
        mod_compound_peaks = [(unknownCompound.spectrum.mz_key[i], unknownCompound.spectrum.intensity[i]) for i in range(len(unknownCompound.spectrum.mz_key))]
        matched_peaks = self.get_edge_detail(knownCompound.id, unknownCompound.id)
        if matched_peaks is None:
            matched_peaks = []
        else:
            matched_peaks = matched_peaks.get_matches_pairs()
        
        peaksObj = {
            "main_compound_peaks": main_compound_peaks,
            "mod_compound_peaks": mod_compound_peaks,
            "matched_peaks": matched_peaks,
            "main_precursor_mz": knownCompound.spectrum.precursor_mz,
            "mod_precursor_mz": unknownCompound.spectrum.precursor_mz,
        }

        fragmentsObj = {
            "frags_map": knownCompound.spectrum.peak_fragment_dict,
            "structure": knownCompound.structure,
            "peaks": main_compound_peaks,
            "Precursor_MZ": knownCompound.spectrum.precursor_mz,
        }

        return peaksObj, fragmentsObj


    def add_adjusted_edge(self, u, v, edgeDetail: EdgeDetail = None, **kwargs):
        """
        Add an edge between two compounds.

        The method will add an edge between two compounds. If the edgeDetail is not passed, the method will align
        the spectra of the compounds using the alignment engine. If the edgeDetail is passed,
        It has to be from the smaller compound to the larger compound.
        
        Parameters
        ----------
        u : str
            The id of the first compound.
        
        v : str
            The id of the second compound.
        
        edgeDetail : EdgeDetail
            The edge detail between the compounds. If not passed, the method will align the spectra of the compounds.
        
        kwargs : dict
            Additional parameters to pass to the alignment engine.
        """
        
        if not self.network.has_node(u):
            raise ValueError(f"{u} is not in the network")
        
        if not self.network.has_node(v):
            raise ValueError(f"{v} is not in the network")
        
        smaller = u if self.network.nodes[u]["compound"].spectrum.precursor_mz <= self.network.nodes[v]["compound"].spectrum.precursor_mz else v
        larger = u if self.network.nodes[u]["compound"].spectrum.precursor_mz > self.network.nodes[v]["compound"].spectrum.precursor_mz else v
        if edgeDetail is None:

            edgeDetail = self.alignmentEngine.align_single(
                self.network.nodes[smaller]["compound"].spectrum,
                self.network.nodes[larger]["compound"].spectrum,
                self.network.nodes[smaller]["compound"].id,
                self.network.nodes[larger]["compound"].id,
                **kwargs)
            
        self.update_edge(smaller, larger, edgeDetail, **kwargs)
    
    def re_align(self, alignmentEngine: AlignmentEngine, **kwargs):
        """
        Re-align the network.
        For each edge in the network, the method will re-align the edges using the alignment engine.
        
        Parameters
        ----------
        alignmentEngine : AlignmentEngine
            The alignment engine to use to align the unknown compound with the known compounds in the network.
            
        kwargs : dict
            Additional parameters to pass to the alignment engine.
        """
        
        if alignmentEngine is None:
            raise ValueError("Alignment engine is required to re-align the network")
        ppm_tolerance = kwargs.get("ppm_tolerance", self.ppm_tolerance)
        kwargs["ppm_tolerance"] = ppm_tolerance
        alignmentEngine.align(self.network, **kwargs)
    
    def re_annotate(self, annotationEngine: AnnotationEngine, **kwargs):
        """
        Annotate the network.
        For each node (Compound) in the network, the method will annotate the node using the annotation engine.
        
        Parameters
        ----------
        annotationEngine : AnnotationEngine
            The annotation engine to use to annotate the compounds in the network.
        
        kwargs : dict
            Additional parameters to pass to the annotation engine.
        """
        
        if annotationEngine is None:
            raise ValueError("Annotation engine is required to annotate the network")
        annotationEngine.annotate(self.network, annotate_all = True, **kwargs)

    def solve(self, unknown: str, **kwargs):
        """
        Solve the network.

        The method will solve the network for an unknown compound
        
        Parameters
        ----------
        unknown : str
            The id of the unknown compound in the network.
        
        alignmentEngine : AlignmentEngine
            The alignment engine to use to align the unknown compound with the known compounds in the network.
            
        kwargs : dict
            Additional parameters.
            
        
        Returns
        -------
        dict
            A dictionary with the following keys:
            - probabilities: the probabilities of the atoms in the unknown compound.
        
        """

        if self.network.in_degree(unknown) + self.network.out_degree(unknown) == 0:
            raise ModiFinderNotSolvableError(
                f"{unknown} is not solveable due to lack of known compounds"
            )

        if self.network.in_degree(unknown) + self.network.out_degree(unknown) > 1:
            # TODO: implement the logic to solve the network when there are multiple known compounds
            raise ModiFinderNotImplementedError(
                f"{unknown} is not solveable due to multiple known compounds"
            )

        if self.network.in_degree(unknown) == 1:
            # TODO: implement the logic to solve the network when there is only one known compound
            pass


    def update_edge(self, u, v, edgeDetail: EdgeDetail = None, **kwargs):
        """
        Update the edge between two compounds.

        The method will update the edge between two compounds.
        """
        if not self.network.has_edge(u, v):
            self.network.add_edge(u, v, edgeDetail=None)
        
        if edgeDetail:
            self.network[u][v]["edgeDetail"] = edgeDetail
        else:
            if not u.spectrum.precursor_mz <= v.spectrum.precursor_mz:
                raise ValueError(
                    f"the edge between {u.id} and {v.id} must be from the smaller compound to the larger compound"
                )
            else:
                raise NotImplementedError("EdgeDetail must be provided to update the edge")
                self.network[u][v]["edgeDetail"] = EdgeDetail()
                

        # update the key in the edgeDetail if passed in kwargs
        for key in self.network[u][v]["edgeDetail"].__dict__:
            if key in kwargs:
                setattr(self.network[u][v]["edgeDetail"], key, kwargs[key])
                
    
    def add_neighbor(self, compound: Compound, neighbor: str, edgeDetail: EdgeDetail = None, **kwargs):
        """
        Add a neighbor to a compound.
        
        The method will add a neighbor to a compound. If the edgeDetail is not passed, the method will align
        the spectra of the compound and the neighbor using the alignment engine. If the edgeDetail is passed,
        It has to be from the smaller compound to the larger compound.
        
        Parameters
        ----------
        compound : Compound
            The compound to add the neighbor to.
        
        neighbor : str
            The id of the neighbor compound.
        
        edgeDetail : EdgeDetail
            The edge detail between the compound and the neighbor.
        
        kwargs : dict
            Additional parameters to pass to the alignment engine.
        """
        
        if not self.network.has_node(neighbor):
            raise ValueError(f"{neighbor} is not in the network")
        
        # check if compound doesn't exist in the network, add it
        if not self.network.has_node(compound.id):
            self.network.add_node(compound.id, compound=compound)
        
        smaller = compound if compound.spectrum.precursor_mz <= self.network.nodes[neighbor]["compound"].spectrum.precursor_mz else self.network.nodes[neighbor]["compound"]
        larger = compound if compound.spectrum.precursor_mz > self.network.nodes[neighbor]["compound"].spectrum.precursor_mz else self.network.nodes[neighbor]["compound"]
        if edgeDetail is None:
            edgeDetail = self.alignmentEngine.align_single(smaller.spectrum, larger.spectrum, **kwargs)
        
        self.update_edge(smaller.id, larger.id, edgeDetail, **kwargs)
        
    
    def add_unknown(self, unknown: Compound, **kwargs) -> str:
        """
        Add an unknown compound to the network.
        
        The method will add an unknown compound to the network.
        
        Parameters
        ----------
        unknown : Compound
            The unknown compound to add to the network.
        
        kwargs : dict
            Additional parameters to pass for the conversion.
        
        Returns
        -------
        str
            The id of the unknown compound.
        """
        
        unknown = convert.to_compound(data=unknown, is_known = False, **kwargs)
        
        if unknown.is_known:
            raise ValueError("Unknown compound must have is_known attribute set to False")
        
        if not self.network.has_node(unknown.id):
            self.network.add_node(unknown.id, compound=unknown)
        
        if not self.unknowns:
            self.unknowns = [unknown]
        else:
            self.unknowns.append(unknown)
        
        return unknown.id
        
    
    def generate_probabilities(self, known_id = None, unknown_id = None, shifted_only = False, CI = False, CPA = True, CFA = True, CPE = True, **kwargs):
        
        if unknown_id is None:
            unknown_id = self.unknowns[0].id
        
        if known_id is None and unknown_id is not None:
            known_id = self.knowns[0].id
        
        if known_id is None or unknown_id is None:
            raise ValueError("Both known and unknown compound ids must be specified")
        
        edgeDetail = self.get_edge_detail(known_id, unknown_id)
        if edgeDetail is None:
            raise ValueError(f"No edge detail found between {known_id} and {unknown_id}")

        # Determine if 'known_id' was the first or second spectrum in the alignment
        # CosineAlignmentEngine stores these IDs during alignment
        if edgeDetail.start_compound_id == known_id:
            shifted_peaks_in_known = [m.first_peak_mz for m in edgeDetail.matches if m.match_type == MatchType.shifted]
            unshifted_peaks_in_known = [m.first_peak_mz for m in edgeDetail.matches if m.match_type == MatchType.unshifted]
        elif edgeDetail.end_compound_id == known_id:
            shifted_peaks_in_known = [m.second_peak_mz for m in edgeDetail.matches if m.match_type == MatchType.shifted]
            unshifted_peaks_in_known = [m.second_peak_mz for m in edgeDetail.matches if m.match_type == MatchType.unshifted]
        else:
            raise ValueError(f"EdgeDetail does not contain compound ID for {known_id}")

        known_compound = self.network.nodes[known_id]["compound"]
        assert known_compound.is_known, f"Compound {known_id} is not a known compound"
        assert known_compound.id == known_id, f"Compound id {known_compound.id} does not match known_id {known_id}"

        positive_contributions = known_compound.calculate_contributions_by_peak_id(shifted_peaks_in_known, CI=CI, CPA=CPA, CFA=CFA, CPE=CPE)
        
        if not shifted_only:
            negative_contributions = known_compound.calculate_contributions_by_peak_id(unshifted_peaks_in_known, CI=CI, CPA=CPA, CFA=CFA, CPE=CPE)
        else:
            negative_contributions = [0 for i in range(len(known_compound.structure.GetAtoms()))]

        probabilities = np.zeros(len(known_compound.structure.GetAtoms()))
        assert len(positive_contributions) == len(negative_contributions) == len(probabilities), "The length of contributions and probabilities must be the same"
        for i in range(len(positive_contributions)):
            probabilities[i] = positive_contributions[i] - negative_contributions[i]
        
        if np.min(probabilities) < 0:
            probabilities = probabilities - np.min(probabilities)
            # mask only the atoms with positive contribution
            # probabilities = [probabilities[i] if positive_contributions[i] > 0 else 0 for i in range(len(probabilities))]
        if np.sum(probabilities) > 0:
            probabilities = probabilities / np.sum(probabilities)

        probabilities = general_utils.power_prob(probabilities)

        return probabilities
    
    def get_edge_detail(self, id1, id2) -> EdgeDetail:
        """Returns the edge detail between two compounds. If the backward edge exists, the matching data is replaced

        Parameters
        ----------
        id1 : str
            The id of the first compound.
        
        id2 : str
            The id of the second compound.
        
        Returns
        -------
        EdgeDetail
            The edge detail between the two compounds.
        
        Raises
        ------
        ValueError
            If the compounds are not connected in the network.
        """
        
        if self.network.has_edge(id1, id2):
            return self.network[id1][id2]["edgeDetail"]
        elif self.network.has_edge(id2, id1):
            edge_detail = self.network[id2][id1]["edgeDetail"]
            new_edge_detail = edge_detail.copy()
            new_edge_detail.reverse_match()
            return new_edge_detail
        else:
            raise ValueError(f"Compounds {id1} and {id2} are not connected in the network")
        
        
    
    
    def draw_prediction(self, probabilities, known_id, **kwargs):
        
        return mf_vis.draw_molecule_heatmap(self.network.nodes[known_id]["compound"].structure, probabilities, **kwargs)
    
    def draw_alignment(self, id1, id2, **kwargs):
        """
        Draw the alignment between two compounds.
        
        The method will draw the alignment between two compounds. The compounds must be in the network and connected by an edge.
        
        See Also
        --------
        visualizer.draw_alignment
        
        Parameters
        ----------
        id1 : str
            The id of the first compound.
        
        id2 : str
            The id of the second compound.
        
        kwargs : dict
            Additional parameters to pass to the visualizer.
        """
        
        if not self.network.has_edge(id1, id2):
            if not self.network.has_edge(id2, id1):
                raise ValueError(f"Compounds {id1} and {id2} are not connected in the network")
            else:
                edgeDetail = self.network[id2][id1]["edgeDetail"]
                edgeDetailCopy = edgeDetail.copy()
                edgeDetailCopy.reverse_match()
                if edgeDetail is None:
                    matched_peaks = []
                else:
                    matched_peaks = edgeDetailCopy.get_matches_pairs()
        else:
            edgeDetail = self.network[id1][id2]["edgeDetail"]
            if edgeDetail is None:
                matched_peaks = []
            else:
                matched_peaks = edgeDetail.get_matches_pairs()

        s1 = self.network.nodes[id1]["compound"].spectrum
        s2 = self.network.nodes[id2]["compound"].spectrum
        
        return mf_vis.draw_alignment([s1, s2], [matched_peaks], **kwargs)
        

    def _get_unknown(self):
        if len(self.unknowns) > 1:
            raise ValueError("More than one unknown compound found in the network. Please specify the unknown compound id")
        unknown_id = self.unknowns[0].id
        return unknown_id
    
    def _get_known_neighbor(self, unknown_id):
        neighbors = list(self.network.predecessors(unknown_id)) + list(self.network.successors(unknown_id))
        if len(neighbors) > 1:
            raise ValueError("More than one known compound found in the network. Please specify the known compound id")
        known_id = neighbors[0]
        
        return known_id
    
    def get_neighbors(self, node_id):
        neighbors = list(self.network.predecessors(node_id)) + list(self.network.successors(node_id))
        return neighbors
    
    def get_meta_data(self, known_id, unknown_id):
        known_compound = self.network.nodes[known_id]["compound"]
        unknown_compound = self.network.nodes[unknown_id]["compound"]
        
        known_meta = known_compound.get_meta_data()
        unknown_meta = unknown_compound.get_meta_data()
        
        result = {}
        for key in known_meta:
            result[key+"_main"] = known_meta[key]
        
        for key in unknown_meta:
            result[key+"_modified"] = unknown_meta[key]
        
        result["delta_mass"] = unknown_compound.spectrum.precursor_mz - known_compound.spectrum.precursor_mz
        result["is_addition"] = 1 if result["delta_mass"] > 0 else -1
        result["num_helpers"] = len(list(self.network.predecessors(known_id)) + list(self.network.successors(known_id))) - 1

        
        try:
            edgeDetail = self.get_edge_detail(known_id, unknown_id)
            result.update(edgeDetail.get_meta_data())
            shifted_peaks = edgeDetail.get_single_type_matches(MatchType.shifted)
            unshifted_peaks = edgeDetail.get_single_type_matches(MatchType.unshifted)
            
            shifted_peaks = [peak[0] for peak in shifted_peaks]
            unshifted_peaks = [peak[0] for peak in unshifted_peaks]
            
            shifted_annotated_ratio, shifted_annotated_ambiguity = known_compound.calculate_peak_annotation_ambiguity(shifted_peaks)
            unshifted_annotated_ratio, unshifted_annotated_ambiguity = known_compound.calculate_peak_annotation_ambiguity(unshifted_peaks)
            
            result["shifted_annotated_ratio"] = shifted_annotated_ratio
            result["shifted_annotated_ambiguity"] = shifted_annotated_ambiguity
            result["unshifted_annotated_ratio"] = unshifted_annotated_ratio
            result["unshifted_annotated_ambiguity"] = unshifted_annotated_ambiguity
            
        except Exception as e:
            print(f"Error getting edge detail: {e}")
            pass
        
        return result
