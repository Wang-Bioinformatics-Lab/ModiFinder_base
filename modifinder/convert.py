import warnings
import modifinder as mf
import modifinder.utilities.network as network
from modifinder.utilities.general_utils import parse_data_to_universal
import json
from copy import deepcopy

def to_compound(data = None, use_object=None, **kwargs):
    """Make a Compound object from the data
    
    Parameters
    ----------
    data: object to be converted

        Current supported types are:
         Compound object (return the same object, for copying you can pass use_object or use .copy() method)
         USI string
         dictionary-of-data
    
    use_object: object, optional
        If a Compound object is passed, this object will be used to create the new object.
    
    kwargs: keyword arguments
        If no data is passed, the keyword arguments will be used to create the object.
    """
    if data is None:
        data = parse_data_to_universal(kwargs)
        try:
            if use_object:
                compound = use_object
                compound.clear()
                compound.update(**data)
            else:
                compound = mf.Compound()
                compound.update(**data)
            return compound
        except Exception as err:
            raise mf.ModiFinderError("Input data is not a valid dictionary.") from err

    # Compound Object
    if hasattr(data, "spectrum"):
        try:
            if use_object:
                compound = use_object
                compound.clear()
                data_dict = compound_to_dict(data)
                data_dict.update(kwargs)
                compound.update(data_dict)
            else:
                compound = data
            return compound

        except Exception as err:
            raise mf.ModiFinderError("Input data is not a valid Compound object.") from err
    
    # USI
    if isinstance(data, str):
        print("in convert, data is a string", data)
        try:
            
            data = network.get_data(data)
            data.update(kwargs)
            if use_object:
                compound = use_object
                compound.clear()
                compound.update(**data)
            else:
                compound = mf.Compound()
                compound.update(**data)
            return compound

        except Exception as err:
            raise mf.ModiFinderError("Input data is not a valid USI string." + str(err)) from err
    
    # Dictionary
    if isinstance(data, dict):
        data = parse_data_to_universal(data)
        data.update(kwargs)
        try:
            if use_object:
                compound = use_object
                compound.clear()
                compound.update(**data)
            else:
                compound = mf.Compound(**data)
            return compound
        except Exception as err:
            raise mf.ModiFinderError("Input data is not a valid dictionary.") from err
        

def compound_to_dict(compound):
    """Convert a Compound object to a dictionary"""
    return compound.__dict__


def to_spectrum(data = None, use_object=None, needs_parse = True, **kwargs):
    """Make a Spectrum object from the data
    
    Parameters
    ----------
    data: object to be converted

        Current supported types are:
         Spectrum object (return the same object, for copying you can pass use_object or use .copy() method)
         USI string
         dictionary-of-data
    
    use_object: object, optional
        If a Spectrum object is passed, this object will be used to create the new object.
    
    needs_parse: bool, default is True
        If True, the dict data will be parsed to a universal format
    """
    if data is None:
        if needs_parse:
            data = parse_data_to_universal(kwargs)
        else:
            data = kwargs
        try:
            if use_object:
                spectrum = use_object
                spectrum.clear()
                spectrum.update(**data)
            else:
                spectrum = mf.Spectrum(incoming_data=data)
            return spectrum
        except Exception as err:
            raise mf.ModiFinderError("Input data is not a valid dictionary.") from err

    # Spectrum Object
    if hasattr(data, "mz"):
        if needs_parse:
            additional_data = parse_data_to_universal(kwargs)
        else:
            additional_data = kwargs
        try:
            if use_object:
                spectrum = use_object
                spectrum.clear()
                spectrun_dict = spectrum_to_dict(data)
                spectrun_dict.update(additional_data)
                spectrum.update(**spectrun_dict)
            else:
                spectrum = data
            return spectrum

        except Exception as err:
            raise mf.ModiFinderError("Input data is not a valid Spectrum object.") from err
    
    # USI
    if isinstance(data, str):
        try:
            data = network.get_data(data)
            data = parse_data_to_universal(data)
            data.update(kwargs)
            data = parse_data_to_universal(data)
            if use_object:
                spectrum = use_object
                spectrum.clear()
                spectrum.update(**data)
            else:
                spectrum = mf.Spectrum(**data)
            return spectrum

        except Exception as err:
            raise mf.ModiFinderError("Input data is not a valid USI string.") from err
    
    # Dictionary
    if isinstance(data, dict):
        if needs_parse:
            data = parse_data_to_universal(data)
        try:
            # add kwargs to data
            data.update(kwargs)
            if use_object:
                spectrum = use_object
                spectrum.clear()
                spectrum.update(**data)
            else:
                spectrum = mf.Spectrum(**data)
            return spectrum
        except Exception as err:
            raise mf.ModiFinderError("Input data is not a valid dictionary.") from err
    
    if isinstance(data, list):
        new_data = dict()
        new_data["mz"] = [x[0] for x in data]
        new_data["intensity"] = [x[1] for x in data]
        try:
            if use_object:
                spectrum = use_object
                spectrum.clear()
                spectrum.update(**new_data)
            else:
                spectrum = mf.Spectrum(**new_data)
            return spectrum
        except Exception as err:
            raise mf.ModiFinderError("Input data is not a valid list.") from err
    
    raise mf.ModiFinderError("Input data is not a valid object.")
    

def spectrum_to_dict(spectrum):
    """Convert a Spectrum object to a dictionary"""
    return spectrum.__dict__