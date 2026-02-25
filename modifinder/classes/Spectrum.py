import json
from modifinder.utilities.gnps_types import adduct_mapping
import modifinder.utilities.general_utils as general_utils
from .. import convert
import numpy as np
import bisect
import uuid

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
    adduct_mass: float
        The adduct mass.
    ms_level: int
        The ms level, default is 2.
    instrument: str, optional
        The instrument used.
    ms_mass_analyzer: str, optional
        The mass analyzer used.
    ms_dissociation_method: str, optional
        The dissociation method used.
    spectrum_id: str, optional
        The spectrum id. if not provided, it will be generated.
    peak_fragment_dict: dict, optional
        A dictionary mapping peaks to fragments
    
    Examples
    --------
    
    """
    
    def __init__(self, incoming_data=None, **kwargs):
        """Constructor for the Spectrum class.

        the spectrum class can be initialized with a dictionary of data or with the individual values.

        Parameters
        ----------
            incoming_data : Input data (optional, default is None).
                The data to initialize the Spectrum object. The data can be a dictionary, a usi, or a Spectrum object.
        """
        self.mz_key = None   # TODO: Finish implementing tomorrow, this is the mz values multiplied by 1e6 and converted to integers for accurate ID
        self.intensity = None
        self.precursor_mz = None
        self.precursor_charge = None
        self._adduct = None
        self._adduct_mass = None
        self.ms_level = None
        self.instrument = None
        self.ms_mass_analyzer = None
        self.ms_dissociation_method = None
        self.spectrum_id = None
        self.peak_fragment_dict = {}
        
        if incoming_data is None and len(kwargs) == 0:
            return
        
        if incoming_data is not None:
            kwargs['normalize_peaks'] = True
            convert.to_spectrum(incoming_data, self, **kwargs)
        else:
            self.update(**kwargs)

    @property
    def adduct(self):
        return self._adduct
    
    @adduct.setter
    def adduct(self, value):
        self._adduct = adduct_mapping.get(value, value)
        if self._adduct is not None:
            try:
                self._adduct_mass = general_utils.get_adduct_mass(self._adduct)
            except Exception as e:
                if not self.ignore_adduct_format:
                    raise e
                else:
                    self._adduct_mass = None
        else:
            self._adduct_mass = 0
    
    @property
    def adduct_mass(self):
        return self._adduct_mass

    def update(self, peaks = None, peaks_json = None, mz=None, intensity=None, precursor_mz=None, precursor_charge=None, 
               _adduct = None, adduct=None, ms_level=None, instrument=None, ms_mass_analyzer=None, 
               ms_dissociation_method=None, spectrum_id=None, peak_fragment_dict: dict = None, ignore_adduct_format = False, **kwargs):
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
            keep_top_k (int): If not None, only keep the top k peaks.
            peak_fragment_dict (dict): A dictionary mapping peaks to fragments
            ignore_adduct_format (bool): If True, if the adduct format is not recognized, it will not throw an error.
            **kwargs: Additional keyword arguments.
        """
        self.ignore_adduct_format = ignore_adduct_format
        if peaks_json is not None:
            peaks = json.loads(peaks_json)
        if peaks is not None:
            self.mz_key = [int(peak[0] * 1e6) for peak in peaks]
            self.intensity = [peak[1] for peak in peaks]
        self.mz_key = np.array([int(x * 1e6) for x in mz]) if mz is not None else self.mz_key
        self.intensity = np.array(intensity) if intensity is not None else self.intensity
        self.precursor_mz = float(precursor_mz) if precursor_mz is not None else self.precursor_mz
        self.precursor_charge = int(float(precursor_charge)) if precursor_charge is not None else self.precursor_charge
        self.adduct = _adduct if _adduct is not None else self.adduct
        self.adduct = adduct if adduct is not None else self.adduct
        self.ms_level = int(ms_level) if ms_level is not None else self.ms_level
        self.instrument = instrument if instrument is not None else self.instrument
        self.ms_mass_analyzer = ms_mass_analyzer if ms_mass_analyzer is not None else self.ms_mass_analyzer
        self.ms_dissociation_method = ms_dissociation_method if ms_dissociation_method is not None else self.ms_dissociation_method
        self.spectrum_id = spectrum_id if spectrum_id is not None else self.spectrum_id
        self.peak_fragment_dict = peak_fragment_dict if peak_fragment_dict is not None else self.peak_fragment_dict
        
        if self.mz_key is not None:
            self.mz_key, self.intensity = zip(*sorted(zip(self.mz_key, self.intensity)))
        
        self.normalize_peaks()
        
        if self.spectrum_id is None:
            self.spectrum_id = str(uuid.uuid4())

        # make sure types are correct
        if self.precursor_mz is not None:
            self.precursor_mz = float(self.precursor_mz)
        if self.precursor_charge is not None:
            self.precursor_charge = int(self.precursor_charge)
            if self.precursor_charge == 0:
                self.precursor_charge = 1


    def __str__(self):
        object_dict = self.__dict__
        to_delete = [keys for keys in object_dict.keys() if object_dict[keys] is None]
        for key in to_delete:
            del object_dict[key]
        return json.dumps(object_dict, indent=4)
    
    
    def clear(self):
        """Clear the Spectrum object."""
        self.mz_key = None
        self.intensity = None
        self.precursor_mz = None
        self.precursor_charge = None
        self._adduct = None
        self._adduct_mass = None
        self.ms_level = None
        self.instrument = None
        self.ms_mass_analyzer = None
        self.ms_dissociation_method = None
        self.spectrum_id = None
        self.peak_fragment_dict = {}
    
    
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
        mzs = self.mz_key / 1e6
        
        min_range = max(mz-mz_tolerance, mz - (mz * ppm_tolerance / 1e6))
        max_range = min(mz+mz_tolerance, mz + (mz * ppm_tolerance / 1e6))
        
        # Find the leftmost index where min_val could be inserted
        left_index = bisect.bisect_left(mzs, min_range)
        # Find the rightmost index where max_val could be inserted
        right_index = bisect.bisect_right(mzs, max_range)
        # Return the range of indices between left_index and right_index (exclusive of right_index)
        return list(range(left_index, right_index))
    
    
