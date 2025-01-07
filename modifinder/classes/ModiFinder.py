"""Base class for ModiFinder.

This class is used to create a ModiFinder object.
The object can be used to get information about unknown compounds in the network using the known compounds.
"""

from modifinder import convert as convert
from modifinder.classes.Compound import Compound
from modifinder.classes.EdgeDetail import EdgeDetail, MatchType
from modifinder.engines.Abtracts import AlignmentEngine, AnnotationEngine
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
        and the "edgedetail" attribute of an edge is an EdgeDetail object.

    unknowns : list
        A list of compound ids that are unknown in the network.

    See Also
    --------
    Compound
    EdgeDetail
    Engines

    Examples
    --------
    """

    def __init__(
        self,
        knownCompond: Compound = None,
        unknownCompound: Compound = None,
        edgeDetail: EdgeDetail = None,
        helpers: list = [],
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
        knownCompond : known compound (not optional in Use Case 1, ignored in Use Case 2)
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
            and the "edgedetail" attribute of an edge is an *EdgeDetail* object.

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
        Engines

        Examples
        --------

        """

        self.network = None
        self.unknonws = None
        self.ppm_tolerance = ppm_tolerance
        self.args = kwargs
        
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
            knownCompond = convert.to_compound(data=knownCompond, **kwargs)
            unknownCompound = convert.to_compound(data=unknownCompound, **kwargs)
            unknownCompound.is_known = False
            self.network.add_node(knownCompond.id, compound=knownCompond)
            self.network.add_node(unknownCompound.id, compound=unknownCompound)

            self.add_adjusted_edge(knownCompond.id, unknownCompound.id, edgeDetail, **kwargs)

            self.unknowns = [unknownCompound.id]
        
        if helpers is not None:
            for helper in helpers:
                try:
                    helper = convert.to_compound(data=helper)
                    self.network.add_node(helper.id, compound=helper)
                    self.add_adjusted_edge(helper.id, knownCompond.id)
                except Exception as e:
                    print(f"Error adding helper compound: {e}")
                    raise e
        
        if should_align:
            self.re_align(self.alignmentEngine, **kwargs)
        
        if should_annotate:
            self.re_annotate(self.annotationEngine, **kwargs)
    
    
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
            The edge detail between the compounds.
        
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
            edgeDetail = self.alignmentEngine.single_align(self.network.nodes[smaller]["compound"].spectrum, self.network.nodes[larger]["compound"].spectrum, **kwargs)
        
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
            self.network[u][v]["edgedetail"] = edgeDetail
        else:
            if not u.spectrum.precursor_mz <= v.spectrum.precursor_mz:
                raise ValueError(
                    f"the edge between {u.id} and {v.id} must be from the smaller compound to the larger compound"
                )
            else:
                self.network[u][v]["edgedetail"] = EdgeDetail()
                

        # update the key in the edgedetail if passed in kwargs
        for key in self.network[u][v]["edgedetail"].__dict__:
            if key in kwargs:
                setattr(self.network[u][v]["edgedetail"], key, kwargs[key])
                
    
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
            edgeDetail = self.alignmentEngine.single_align(smaller.spectrum, larger.spectrum, **kwargs)
        
        self.update_edge(smaller.id, larger.id, edgeDetail, **kwargs)
        
    
    def generate_probabilities(self, known_id = None, unknown_id = None, shifted_only = False, CI = False, CPA = True, CFA = True, CPE = True):
        
        if unknown_id is None:
            unknown_id = self._get_unknown()
        
        if known_id is None and unknown_id is not None:
            known_id = self._get_known_neighbor(unknown_id)
        
        if known_id is None or unknown_id is None:
            raise ValueError("Both known and unknown compound ids must be specified")
        
        # check if unknown is connected to known
        if self.network.has_edge(known_id, unknown_id):
            edgeDetail = self.network[known_id][unknown_id]["edgedetail"]
            if edgeDetail is None:
                raise ValueError(f"Edge between {known_id} and {unknown_id} does not have edge details")
        
            shifted_peaks_in_known = [match.first_peak_index for match in edgeDetail.matches if match.match_type == MatchType.shifted]
            unshifted_peaks_in_known = [match.first_peak_index for match in edgeDetail.matches if match.match_type == MatchType.unshifted]
        elif self.network.has_edge(unknown_id, known_id):
            edgeDetail = self.network[unknown_id][known_id]["edgedetail"]
            if edgeDetail is None:
                raise ValueError("No edge detail found between the known and unknown compounds")    
            shifted_peaks_in_known = [match.second_peak_index for match in edgeDetail.matches if match.match_type == MatchType.shifted]
            unshifted_peaks_in_known = [match.second_peak_index for match in edgeDetail.matches if match.match_type == MatchType.unshifted]
        else:
            raise ValueError("No edge found between the known and unknown compounds")
        known_compound = self.network.nodes[known_id]["compound"]
        positive_contributions = known_compound.calculate_contributions(shifted_peaks_in_known, CI=CI, CPA=CPA, CFA=CFA, CPE=CPE)
        if not shifted_only:
            negative_contributions = known_compound.calculate_contributions(unshifted_peaks_in_known, CI=CI, CPA=CPA, CFA=CFA, CPE=CPE)
        else:
            negative_contributions = [0 for i in range(len(known_compound.structure.GetAtoms()))]
        
        probabilities = np.zeros(len(known_compound.structure.GetAtoms()))
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
    
    
    def get_edge_detail(self, id1, id2):
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
            return self.network[id1][id2]["edgedetail"]
        elif self.network.has_edge(id2, id1):
            edge_detail = self.network[id2][id1]["edgedetail"]
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
                edgeDetail = self.network[id2][id1]["edgedetail"]
                smallerSpectrum = self.network.nodes[id2]["compound"].spectrum
                largerSpectrum = self.network.nodes[id1]["compound"].spectrum
                if edgeDetail is None:
                    matched_peaks = []
                else:
                    matched_peaks = [(match.second_peak_index, match.first_peak_index) for match in edgeDetail.matches]
        else:
            edgeDetail = self.network[id1][id2]["edgedetail"]
            smallerSpectrum = self.network.nodes[id1]["compound"].spectrum
            largerSpectrum = self.network.nodes[id2]["compound"].spectrum
            if edgeDetail is None:
                matched_peaks = []
            else:
                matched_peaks = [(match.first_peak_index, match.second_peak_index) for match in edgeDetail.matches]
        
        return mf_vis.draw_alignment([smallerSpectrum, largerSpectrum], [matched_peaks], **kwargs)
        
        
        

    def _get_unknown(self):
        if len(self.unknowns) > 1:
            raise ValueError("More than one unknown compound found in the network. Please specify the unknown compound id")
        unknown_id = self.unknowns[0]
        return unknown_id
    
    def _get_known_neighbor(self, unknown_id):
        neighbors = list(self.network.predecessors(unknown_id)) + list(self.network.successors(unknown_id))
        if len(neighbors) > 1:
            raise ValueError("More than one known compound found in the network. Please specify the known compound id")
        known_id = neighbors[0]
        
        return known_id