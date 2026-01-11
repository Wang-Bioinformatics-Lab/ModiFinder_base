from modifinder.utilities import general_utils as mf_gu
from modifinder.classes.Spectrum import Spectrum
import tqdm

def _cluster_spectrums(spectrums, ppm_tolerance = 10, mz_tolerance = 0.1, verbose = False):
    '''
    This function is used to cluster the peaks of the spectrums within a given tolerance
    
    Parameters
    ----------
        spectrums: list of Spectrum objects
            The spectrums to be clustered
        ppm_tolerance: float
            The tolerance in ppm
        mz_tolerance: float
            The tolerance in mz
    
    
    Returns:
    clusters: list of clusters of the mz values, each cluster is a list of [mz, intensity, origin_spectrum_index]
    '''
    
    all_points = []
    for index, spectrum in tqdm.tqdm(enumerate(spectrums), disable = not verbose):
        for mz, intensity in zip(spectrum.mz, spectrum.intensity):
            all_points.append((mz, intensity, index))
    
    # print("all points", len(all_points))

    # 1. Sort by mz
    all_points.sort(key=lambda x: x[0])

    # 2. Sweep
    clusters = []
    if not all_points:
        return clusters

    # Start the first cluster
    first_mz, first_intensity, first_idx = all_points[0]
    clusters.append([[first_mz, first_intensity, first_idx]])
    ref_mz = first_mz

    for (mz, intensity, idx) in all_points[1:]:
        if not mf_gu.is_shifted(mz, ref_mz, ppm_tolerance, mz_tolerance):
            # Add to the last cluster
            clusters[-1].append([mz, intensity, idx])
        else:
            # Create a new cluster
            clusters.append([[mz, intensity, idx]])
            ref_mz = mz

    return clusters


def aggregate_spectrums(spectrums, ppm_tolerance = 10, mz_tolerance = 0.1, consensus_majority_ratio = 1, **kwargs):
    '''
    This function is to create a consensus spectrum from a list of spectrums of the same compound
    all the other components of the generated spectrum such as the charge, ... are taken from the first spectrum.
    
    if the spectrums have different adducts, the adduct mass is removed from the spectrums before aggregation. The final result is reported with no adduct mass.
    Parameters
    ----------
        spectrums: list of Spectrum objects
            The spectrums to be aggregated
        ppm_tolerance: float
            The tolerance in ppm
        mz_tolerance: float
            The tolerance in mz
        consensus_majority_ratio: float
            The ratio of the spectrums that should have the peak to be considered in the consensus
            if is 1, all the spectrums should have the peak to be considered in the consensus
            if is 0, only one spectrum should have the peak to be considered in the consensus
    
    Returns:
    result: Spectrum object
        The aggregated spectrum
    '''
    
    adducts = set([s.adduct for s in spectrums])
    spectrums = [remove_adduct(s) for s in spectrums]
    
    clusters = _cluster_spectrums(spectrums, ppm_tolerance, mz_tolerance)
    # print("clustering done", len(clusters))
    
    aggregate_mz = []
    aggregate_intensity = []
    for cluster in clusters:
        sum_intensity = 0
        sum_mz = 0
        consensus = set()
        for mz, intensity, index in cluster:
            sum_intensity += intensity
            sum_mz += mz
            consensus.add(index)
        
        if len(consensus) < consensus_majority_ratio * len(spectrums):
            continue
        aggregate_mz.append(sum_mz/len(cluster))
        aggregate_intensity.append(sum_intensity/len(cluster))
    
    result = Spectrum(spectrums[0].copy(), mz = aggregate_mz, intensity = aggregate_intensity, **kwargs)
    if len(adducts) == 1:
        result = add_adduct(result, adducts.pop())
    
    return result


def refine_consensus(spectrum, spectrums, ppm_tolerance = 10, mz_tolerance = 0.1, consensus_majority_ratio = 0):
    '''
    This function is used to refine a spectrum based on the consensus of the other spectrums
    any peak that is not in the consensus of the other spectrums is removed.
    
    If the spectrums are with different adducts, the spectrums are adjusted to remove the adduct mass. The final spectrum is adjusted to have the correct adduct mass.
    
    Parameters
    ----------
        spectrum: Spectrum object
            The spectrum to be refined
        spectrums: list of Spectrum objects
            The spectrums to be used to refine the spectrum
        ppm_tolerance: float
            The tolerance in ppm
        mz_tolerance: float
            The tolerance in mz
        consensus_majority_ratio: float (0-1) [default: 0]
            The ratio of the spectrums that should have the peak in order to keep it, 0 means none of the spectrums need to have it
            1 means all of the spectrums need to have it.
            
    Returns:
        A tuple of two elements:
    result: Spectrum object
        The refined spectrum
    appearances: tuple of lists
        The first list contains the number of spectrums that have the peak after refining
        The second list contains the number of spectrums that have the peak before refining
    '''
    # Create an empty list to store the refined consensus
    
    refined_mz = []
    refined_intensity = []
    new_appearances = []
    old_appearances = []
    
    adduct = spectrum.adduct
    spectrum = remove_adduct(spectrum)
    spectrums = [remove_adduct(s) for s in spectrums]
    
    for mz, intensity in zip(spectrum.mz, spectrum.intensity):
        found_count = 0
        for other_spectrum in spectrums:
            peak_indexes = other_spectrum.get_peak_indexes(mz, ppm_tolerance = ppm_tolerance, mz_tolerance = mz_tolerance)
            if len(peak_indexes) > 0:
                found_count += 1
        
        old_appearances.append(found_count)
        
        if found_count >= consensus_majority_ratio * len(spectrums):
            refined_mz.append(mz)
            refined_intensity.append(intensity)
            new_appearances.append(found_count)
    
    result = Spectrum(spectrum.copy(), mz = refined_mz, intensity = refined_intensity)
    result = add_adduct(result, adduct)
    
    return result, (new_appearances, old_appearances)


def remove_adduct(spectrum):
    '''
    This function is used to remove the adduct mass from the spectrum
    
    Parameters
    ----------
        spectrum: Spectrum object
            The spectrum to remove the adduct mass from its mz values
    
    Returns:
    result: Spectrum object
        The spectrum without the adduct mass
    '''
    
    result = spectrum.copy()
    adduct_mass = spectrum.adduct_mass
    result.mz = [mz - adduct_mass for mz in result.mz]
    if result.precursor_mz is not None:
        result.precursor_mz -= adduct_mass
    result.adduct = None
    
    return result


def add_adduct(spectrum, adduct):
    '''
    This function is used to add the adduct mass to the spectrum
    
    Parameters
    ----------
        spectrum: Spectrum object
            The spectrum to add the adduct mass to
        adduct: str
            The adduct to add to the spectrum, should be in the supported adducts list.
    
    Returns:
    result: Spectrum object
        The spectrum with the adduct mass
    '''
    
    result = spectrum.copy()
    adduct_mass = mf_gu.get_adduct_mass(adduct)
    result.mz = [mz + adduct_mass for mz in result.mz]
    if result.precursor_mz is not None:
        result.precursor_mz += adduct_mass
    result.adduct = adduct
    
    return result
