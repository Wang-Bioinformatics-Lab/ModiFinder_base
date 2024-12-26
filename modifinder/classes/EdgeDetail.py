from dataclasses import dataclass, field
from typing import List, Tuple
from enum import Enum
import networkx as nx

class MatchType(Enum):
    shifted = 1
    unshifted = 2

class Match:
    """
    Class for Match
    """
    def __init__(self, first_peak_index:int, second_peak_index:int, match_type: MatchType):
        """
        Initialize the Match object.

        Parameters
        ----------
        first_peak_index : int
            Index of the first peak
        second_peak_index : int
            Index of the second peak
        match_type : MatchType
            Type of the match

        Customization
        -------------
        if you need matches with more information, create a new class and inherit from this class.
        """
        self.first_peak_index = first_peak_index
        self.second_peak_index = second_peak_index
        self.match_type = match_type
    
    def __repr__(self):
        return f"Match({self.first_peak_index}, {self.second_peak_index}, {self.match_type})"

class EdgeDetail:
    """
    Class for Edge Details
    """
    def __init__(self, number_of_modifications: int = -1, match_score: float = 0, matches: List[Match] = None, start_spectrum_id: str = None, end_spectrum_id: str = None):
        """
        Initialize the EdgeDetail object.

        Parameters
        ----------
        number_of_modifications : int
            Number of modifications, -1 for unknown
        match_score : float
            Match score, how well the two spectra match
        matches : List[Match]
            List of matches, each match is a tuple of two peak indices and the match type. It is 
            important to note that match has directionality. The first peak index is from the first 
            node of the edge and the second peak index is from the second node of the edge.
        """
        self.start_spectrum_id = start_spectrum_id
        self.end_spectrum_id = end_spectrum_id
        self.number_of_modifications = number_of_modifications
        self.match_score = match_score
        self.matches = matches if matches else []
    
    def __str__(self):
        return f"EdgeDetail({self.number_of_modifications}, {self.match_score}, {self.matches})"