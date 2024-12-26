# from __future__ import annotations
from abc import ABC, abstractmethod
from modifinder.classes.Spectrum import Spectrum
from modifinder.classes.Compound import Compound
import networkx as nx
from modifinder.classes.EdgeDetail import EdgeDetail
from typing import List, Tuple


# Base class for alignment engines
class AlignmentEngine(ABC):
    @abstractmethod
    def __init__(self, **kwargs):
        pass
    
    @abstractmethod
    def align(self, network:nx.DiGraph,
        mz_tolerance: float = 0.02,
        ppm_tolerance: float = 100.0,
        align_all: bool = True,
        **kwargs):
        """Aligns the spectra in the network
        
        For each edge in the network, aligns the spectrum of the start node with the spectrum of the end node.
        If the edge has already been aligned, and align_all is False, the edge will not be realigned.

        Parameters
        ----------
            network (nx.DiGraph) : The Compound Graph object to align the spectra in.
            mz_tolerance (float, optional) : The mz tolerance in Da for the fragments. Defaults to 0.02Da. 
            ppm_tolerance (float, optional) : The mz tolerance in ppm for the fragments. Defaults to 100.0ppm.
            align_all (bool, optional) : If True, all edges will be aligned. If False, only the edges that have not been aligned will be aligned. Defaults to False.
        """
        pass
    
    @abstractmethod
    def single_align(self, SpectrumTuple1,
                      SpectrumTuple2, 
                      mz_tolerance: float = 0.02, 
                      ppm_tolerance: float = 100.0, 
                      **kwargs):
        """Aligns two spectra, returns the alignment details as an EdgeDetail object

        Parameters
        ----------
            SpectrumTuple1 (SpectrumTuple) : First spectrum
            SpectrumTuple2 (SpectrumTuple) : Second spectrum
            mz_tolerance (float) : Fragment mz tolerance
            ppm_tolerance (float) : Fragment ppm tolerance
            kwargs : additional arguments
        
        Returns
        -------
            EdgeDetail : the edge detail object
        """
        pass

# Base class for annotation engines
class AnnotationEngine(ABC):
    @abstractmethod
    def __init__(self, **kwargs):
        pass
    
    @abstractmethod
    def annotate(self, network, **kwargs):
        pass

    
    @abstractmethod
    def annotate_single(self, compound,  modify_compound: bool = True, **kwargs) -> List[Tuple[int, List[str]]]:
        """
        provides annotation for the peaks in a single compound

        Parameters:
            :compound (Compound): the compound to be annotated
            :modify_compound (bool): whether to modify the passed compound with the annotations
            :kwargs: additional arguments
        
        Returns:
            :Mapping[int, List[str]]: a dictionary with the indices of the peaks as keys and the list of annotations as values
        """
        pass


    @abstractmethod
    def get_fragment_info(self, Compound: Compound, fragment: int) -> Tuple[List[int], List[Tuple[int, int]], str, str]:
        """
        converts a fragment to a SMILES string

        Parameters:
            :Compound (Compound): the compound
            :fragment (int): the fragment
        
        Returns:
            :tuple(atomlist -> List[int], edge_list -> List[Tuple[int, int]], formula -> str, smiles -> str): a tuple containing the atom list, the edge list, the formula and the SMILES string
        """
        pass

# # Base class for prediction engines
# class PredictionEngine(ABC):
#     @abstractmethod
#     def predict(self, network, **kwargs):
#         pass

#     @abstractmethod
#     def confidence(self, network, prediction, **kwargs):
#         pass