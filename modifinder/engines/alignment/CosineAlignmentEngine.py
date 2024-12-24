"""This module contains the CosineAlignmentEngine class for aligning two spectra using cosine similarity."""

import numpy as np
from typing import List, Tuple
from modifinder.engines.Abtracts import AlignmentEngine
from modifinder.classes.Spectrum import Spectrum
from modifinder.classes.EdgeDetail import EdgeDetail, Match, MatchType
from modifinder.utilities.general_utils import is_shifted
from modifinder.exceptions import ModiFinderError
import networkx as nx


def _cosine_fast(
    spec: Spectrum,
    spec_other: Spectrum,
    mz_tolerance: float = None,
    ppm_tolerance: float = None,
    allow_shift: bool = True,
) -> Tuple[float, List[Tuple[int, int]]]:
    """Approximates alignment with the highest cosine score between two spectra.
    
    if both tolerances are provided, both criteria are used to filter the matches. If only one is provided, only that one is used.

    Args:
        spec (Spectrum): First spectrum
        spec_other (Spectrum): Second spectrum
        mz_tolerance (float): Tolerance in Da for fragment m/z values, if None, it is not used.
        ppm_tolerance (float): Tolerance in ppm for fragment m/z values, if None, it is not used.
        allow_shift (bool): _description_

    Returns:
        Tuple[float, List[Tuple[int, int]]]: _description_
    """
    
    if mz_tolerance is None and ppm_tolerance is None:
        raise ModiFinderError("At least one of mz_tolerance or ppm_tolerance must be provided.")
    
    precursor_charge = max(spec.precursor_charge, 1)
    precursor_mass_diff = (
        spec.precursor_mz - spec_other.precursor_mz
    ) * precursor_charge
    # Only take peak shifts into account if the mass difference is relevant.
    num_shifts = 1
    if allow_shift and abs(precursor_mass_diff) >= mz_tolerance:
        num_shifts += precursor_charge
    other_peak_index = np.zeros(num_shifts, np.uint16)
    mass_diff = np.zeros(num_shifts, np.float32)
    for charge in range(1, num_shifts):
        mass_diff[charge] = precursor_mass_diff / charge

    # Find the matching peaks between both spectra.
    peak_match_scores, peak_match_idx = [], []
    for peak_index, (peak_mz, peak_intensity) in enumerate(
        zip(spec.mz, spec.intensity)
    ):
        # Advance while there is an excessive mass difference.
        for cpi in range(num_shifts):
            while other_peak_index[cpi] < len(spec_other.mz) - 1 and (
                peak_mz - mz_tolerance
                > spec_other.mz[other_peak_index[cpi]] + mass_diff[cpi]
            ):
                other_peak_index[cpi] += 1

        # Match the peaks within the fragment mass window if possible.
        for cpi in range(num_shifts):
            index = 0
            other_peak_i = other_peak_index[cpi] + index
            while (
                other_peak_i < len(spec_other.mz)
                and abs(peak_mz - (spec_other.mz[other_peak_i] + mass_diff[cpi]))
                <= mz_tolerance
            ):
                if abs(peak_mz - (spec_other.mz[other_peak_i] + mass_diff[cpi])) <= (
                    ppm_tolerance * peak_mz / 1e6
                ):
                    peak_match_scores.append(
                        peak_intensity * spec_other.intensity[other_peak_i]
                    )
                    peak_match_idx.append((peak_index, other_peak_i))
                index += 1
                other_peak_i = other_peak_index[cpi] + index

    score, peak_matches = 0.0, []
    if len(peak_match_scores) > 0:
        # Use the most prominent peak matches to compute the score (sort in
        # descending order).
        peak_match_scores_arr = np.asarray(peak_match_scores)
        peak_match_order = np.argsort(peak_match_scores_arr)[::-1]
        peak_match_scores_arr = peak_match_scores_arr[peak_match_order]
        peak_match_idx_arr = np.asarray(peak_match_idx)[peak_match_order]
        peaks_used, other_peaks_used = set(), set()
        for peak_match_score, peak_i, other_peak_i in zip(
            peak_match_scores_arr,
            peak_match_idx_arr[:, 0],
            peak_match_idx_arr[:, 1],
        ):
            if peak_i not in peaks_used and other_peak_i not in other_peaks_used:
                score += peak_match_score
                # Save the matched peaks.
                peak_matches.append((peak_i, other_peak_i))
                # Make sure these peaks are not used anymore.
                peaks_used.add(peak_i)
                other_peaks_used.add(other_peak_i)

    return score, peak_matches


class CosineAlignmentEngine(AlignmentEngine):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        
    def align(
        self,
        network: nx.DiGraph,
        mz_tolerance: float = 0.02,
        ppm_tolerance: float = 100.0,
        align_all: bool = False,
        **kwargs,
    ):
        """
        Aligns the spectra in the network using cosine similarity.
        
        Parameters
        ----------
        network : nx.DiGraph
            The Compound Graph object to align the spectra in.
        mz_tolerance : float, optional
            Fragment mz tolerance, by default 0.02
        ppm_tolerance : float, optional
            Fragment ppm tolerance, by default 100.0
        align_all : bool, default False
            if True, all edges will be aligned, if False, only the edges that have not been aligned will be aligned
        kwargs : dict
            additional arguments
        """
        edges = network.edges(data=True)
        for edge in edges:
            start_compound = network.nodes[edge[0]]["compound"]
            end_compound = network.nodes[edge[1]]["compound"]
            if "edgedetail" not in edge[2] or edge[2]["edgedetail"] is None or align_all:
                edge[2]["edgedetail"] = self.single_align(
                    start_compound.spectrum,
                    end_compound.spectrum,
                    mz_tolerance,
                    ppm_tolerance,
                    **kwargs,
                )
                
        

    def single_align(
        self,
        SpectrumTuple1: Spectrum,
        SpectrumTuple2: Spectrum,
        mz_tolerance: float = 0.02,
        ppm_tolerance: float = 100.0,
        **kwargs,
    ) -> EdgeDetail:
        """
        Aligns two spectra using cosine similarity and returns the cosine score and the matched peaks.

        Parameters:
            SpectrumTuple1 (SpectrumTuple): First spectrum
            SpectrumTuple2 (SpectrumTuple): Second spectrum
            mz_tolerance (float): Fragment mz tolerance
            ppm_tolerance (float): Fragment ppm tolerance
            kwargs: additional arguments

        Returns:
            EdgeDetail: the edge detail object
        """
        cosine, matched_peaks = _cosine_fast(
            SpectrumTuple1,
            SpectrumTuple2,
            mz_tolerance,
            ppm_tolerance,
            True,
        )

        Matches = []
        for match in matched_peaks:
            if is_shifted(
                SpectrumTuple1.mz[match[0]],
                SpectrumTuple2.mz[match[1]],
                ppm_tolerance,
                mz_tolerance,
            ):
                Matches.append(Match(match[0], match[1], MatchType.shifted))
            else:
                Matches.append(Match(match[0], match[1], MatchType.unshifted))

        return EdgeDetail(
            match_score=cosine, matches=Matches, number_of_modifications=-1
        )
