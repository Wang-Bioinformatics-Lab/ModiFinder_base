import json
from modifinder.utilities.gnps_types import adduct_mapping
from modifinder import convert
import numpy as np

class Spectrum:
    """A class to represent a spectrum.
    
    Parameters
    ----------
    mz: list
        A list of m/z values.
    intensity: list
        A list of intensity values.
    precursor_mz: float
        The precursor m/z value.
    precursor_charge: int
        The precursor charge.
    adduct: str
        The adduct.
    ms_level: int
        The ms level, default is 2.
    instrument: str, optional
        The instrument used.
    ms_mass_analyzer: str, optional
        The mass analyzer used.
    ms_dissociation_method: str, optional
        The dissociation method used.
    spectrum_id: str, optional
        The spectrum id.
    
    Examples
    --------
    
    """
    
    def __init__(self, incoming_data=None, normalize_peaks = True, **kwargs):
        """Constructor for the Spectrum class.

        the spectrum class can be initialized with a dictionary of data or with the individual values.

        Parameters
        ----------
            incoming_data : Input data (optional, default is None).
                The data to initialize the Spectrum object. The data can be a dictionary, a usi, or a Spectrum object.
            normalize_peaks : bool, default is True.
                If True, the intensity of the peaks will be normalized at initialization.
        """
        self.mz = None
        self.intensity = None
        self.precursor_mz = None
        self.precursor_charge = None
        self.adduct = None
        self.ms_level = None
        self.instrument = None
        self.ms_mass_analyzer = None
        self.ms_dissociation_method = None
        self.spectrum_id = None
        
        if incoming_data is None and len(kwargs) == 0:
            return
        
        if incoming_data is not None:
            convert.to_spectrum(incoming_data, self)

        self.update(normalize_peaks = normalize_peaks, **kwargs)


    def update(self, peaks = None, peaks_json = None, mz=None, intensity=None, precursor_mz=None, precursor_charge=None, 
               adduct=None, ms_level=None, instrument=None, ms_mass_analyzer=None, 
               ms_dissociation_method=None, spectrum_id=None, normalize_peaks = False, ratio_to_base_peak = None, remove_large_peaks = False, **kwargs):
        """Update the Spectrum object with the given values.

        Args:
            peaks (list): A list of peaks in the form of [[mz1, intensity1], [mz2, intensity2], ...].
            peaks_json (str): A json string of peaks.
            mz (list): A list of m/z values.
            intensity (list): A list of intensity values.
            precursor_mz (float): The precursor m/z value.
            precursor_charge (int): The precursor charge.
            adduct (str): The adduct.
            ms_level (int): The ms level.
            instrument (str): The instrument used.
            ms_mass_analyzer (str): The mass analyzer used.
            ms_dissociation_method (str): The dissociation method used.
            spectrum_id (str): The spectrum id.
            normalize_peaks (bool): If True, the intensity of the peaks will be normalized.
            ratio_to_base_peak (float): If None, no filtering is done, if a float number, it removes all the peaks with intensity
            less than ratio times the base peak.
        """
        if peaks_json is not None:
            peaks = json.loads(peaks_json)
        if peaks is not None:
            self.mz = [peak[0] for peak in peaks]
            self.intensity = [peak[1] for peak in peaks]
        self.mz = np.array(mz) if mz is not None else self.mz
        self.intensity = np.array(intensity) if intensity is not None else self.intensity
        self.precursor_mz = float(precursor_mz) if precursor_mz is not None else self.precursor_mz
        self.precursor_charge = int(float(precursor_charge)) if precursor_charge is not None else self.precursor_charge
        self.adduct = adduct_mapping.get(adduct, adduct) if adduct is not None else self.adduct
        self.ms_level = int(ms_level) if ms_level is not None else self.ms_level
        self.instrument = instrument if instrument is not None else self.instrument
        self.ms_mass_analyzer = ms_mass_analyzer if ms_mass_analyzer is not None else self.ms_mass_analyzer
        self.ms_dissociation_method = ms_dissociation_method if ms_dissociation_method is not None else self.ms_dissociation_method
        self.spectrum_id = spectrum_id if spectrum_id is not None else self.spectrum_id
        
        if self.mz is not None:
            self.mz, self.intensity = zip(*sorted(zip(self.mz, self.intensity)))
        
        if normalize_peaks:
            self.normalize_peaks()
        
        if ratio_to_base_peak is not None:
            self.remove_small_peaks(ratio_to_base_peak)
        
        if remove_large_peaks:
            self.remove_larger_than_precursor_peaks()


    def __str__(self):
        object_dict = self.__dict__
        to_delete = [keys for keys in object_dict.keys() if object_dict[keys] is None]
        for key in to_delete:
            del object_dict[key]
        return json.dumps(object_dict, indent=4)
    
    
    def clear(self):
        """Clear the Spectrum object."""
        self.mz = None
        self.intensity = None
        self.precursor_mz = None
        self.precursor_charge = None
        self.adduct = None
        self.ms_level = None
        self.instrument = None
        self.ms_mass_analyzer = None
        self.ms_dissociation_method = None
        self.spectrum_id = None
    
    
    def copy(self):
        """Return a copy of the Spectrum object."""
        copied_spectrum = Spectrum()
        convert.to_spectrum(self, use_object=copied_spectrum, needs_parse=False)
        return copied_spectrum
    
    
    def normalize_peaks(self, change_self = True):
        """l2 Normalize the intensity of the Spectrum object.
        
        Parameters
        ----------
        change_self : bool, default is True
            If True, the intensity of the Spectrum object will be normalized in place.
            If False, a new Spectrum object with the normalized intensity will be returned.
        
        Returns
        -------
        None
            If change_self is True, the intensity of the Spectrum object will be normalized in place.
        Spectrum
            A new Spectrum object with the normalized intensity.
        """
        
        l2_norm = np.linalg.norm(self.intensity)
        new_intensity = [intensity / l2_norm for intensity in self.intensity]
        
        if change_self:
            self.intensity = new_intensity
        else:
            new_spectrum = self.copy()
            new_spectrum.intensity = new_intensity
            return new_spectrum
    
    
    def remove_small_peaks(self, ratio_to_base_peak:float = 0.01, change_self = True):
        """Remove peaks with intensity lower than a given ratio to the base peak.
        
        Parameters
        ----------
        ratio_to_base_peak : float (0, 1), default is 0.01
            The ratio to the base peak.
        change_self : bool, default is True
            If True, the peaks with intensity lower than the given ratio will be removed in place.
            If False, a new Spectrum object with the peaks removed will be returned.
        
        Returns
        -------
        map from old index to new index
            If change_self is True, the peaks with intensity lower than the given ratio will be removed in place and a map from the old index to the new index will be returned.
        (Spectrum, map from old index to new index)
            A new Spectrum object with the peaks removed and a map from the old index to the new index.
        """
        
        base_peak = max(self.intensity)
        new_mz = []
        new_intensity = []
        index_mapping = {}
        for index, intensity in enumerate(self.intensity):
            if intensity >= ratio_to_base_peak * base_peak:
                new_mz.append(self.mz[index])
                new_intensity.append(intensity)
                index_mapping[index] = len(new_mz) - 1
        
        if change_self:
            self.mz = new_mz
            self.intensity = new_intensity
            return index_mapping
        else:
            new_spectrum = self.copy()
            new_spectrum.mz = new_mz
            new_spectrum.intensity = new_intensity
            return new_spectrum, index_mapping
    
    
    def keep_top_k(self, k:int = 100, change_self: bool = True):
        """Keep only the top k peaks in the Spectrum object.
        
        Parameters
        ----------
        k : int, default is 100
            The number of peaks to keep.
        change_self : bool, default is True
            If True, only the top k peaks will be kept in place.
            If False, a new Spectrum object with only the top k peaks will be returned.
        
        Returns
        -------
        map from old index to new index
            If change_self is True, only the top k peaks will be kept in place and a map from the old index to the new index will be returned.
        (Spectrum, map from old index to new index)
            A new Spectrum object with only the top k peaks and a map from the old index to the new index.
        """
        
        top_k_indices = np.argsort(self.intensity)[::-1][:k]
        new_mz = np.array([self.mz[index] for index in top_k_indices])
        new_intensity = np.array([self.intensity[index] for index in top_k_indices])
        # sort the peaks by mz
        new_mz, new_intensity = zip(*sorted(zip(new_mz, new_intensity)))
        index_mapping = {index: new_index for new_index, index in enumerate(top_k_indices)}
        
        if change_self:
            self.mz = new_mz
            self.intensity = new_intensity
            return index_mapping
        else:
            new_spectrum = self.copy()
            new_spectrum.mz = new_mz
            new_spectrum.intensity = new_intensity
            return new_spectrum, index_mapping
        
    def remove_larger_than_precursor_peaks(self, change_self: bool = True):
        """
        Remove peaks that are larger than the precursor m/z value.
        
        Parameters
        ----------
        change_self : bool, default is True
            If True, the peaks that are larger than the precursor m/z value will be removed in place.
            If False, a new Spectrum object with the peaks removed will be returned.
        
        """
        
        new_mz = []
        new_intensity = []
        for mz, intensity in zip(self.mz, self.intensity):
            if mz <= self.precursor_mz:
                new_mz.append(mz)
                new_intensity.append(intensity)
        
        if change_self:
            self.mz = new_mz
            self.intensity = new_intensity
        
        else:
            new_spectrum = self.copy()
            new_spectrum.mz = new_mz
            new_spectrum.intensity = new_intensity
            return new_spectrum
    
    def get_peak_indexes(self, mz, mz_tolerance = 0.02, ppm_tolerance = 40.0, **kwargs):
        """Get the indexes of the peaks within the given m/z tolerance.
        
        Parameters
        ----------
        mz : float
            The m/z value of the peak.
        mz_tolerance : float, optional
            The m/z tolerance for the peak, default is 0.02.
        ppm_tolerance : float, optional
            The ppm tolerance for the peak, default is 40.0.
        
        Returns
        -------
        list
            The indexes of the peaks within the given m/z tolerance.
        """
        
        result = []
        for index, peak_mz in enumerate(self.mz):
            if abs(peak_mz - mz) <= mz_tolerance or abs(peak_mz - mz) <= mz * ppm_tolerance / 1e6:
                result.append(index)    
        return result
    
    

