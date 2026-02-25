import modifinder.utilities.network as network
from modifinder.utilities.general_utils import parse_data_to_universal
from modifinder.exceptions import ModiFinderError
from copy import deepcopy

def _get_compound_class():
    from modifinder.classes.Compound import Compound
    return Compound

def _get_spectrum_class():
    from modifinder.classes.Spectrum import Spectrum
    return Spectrum

def to_compound(data = None, use_object=None, **kwargs):
    """Make a Compound object from the data
    
    Parameters
    ----------
    data: object to be converted

        Current supported types are:
         Compound object (return the same object, for copying you can pass use_object or use .copy() method)
         dictionary-of-data
    
    use_object: object, optional
        If a Compound object is passed, this object will be used to create the new object.
    
    kwargs: keyword arguments
        If no data is passed, the keyword arguments will be used to create the object.
    """
    if data is None:
        data = parse_data_to_universal(kwargs)
        try:
            Compound = _get_compound_class()
            if use_object:
                compound = use_object
                compound.clear()
                compound.update(**data)
            else:
                compound = Compound()
                compound.update(**data)
            return compound
        except Exception as err:
            raise ModiFinderError("Input data is not a valid dictionary. " + str(err)) from err

    # Compound Object
    if isinstance(data, _get_compound_class()):
        return deepcopy(data)
    
    # Dictionary
    if isinstance(data, dict):
        data = parse_data_to_universal(data)
        data.update(kwargs)
        try:
            Compound = _get_compound_class()
            if use_object:
                compound = use_object
                compound.clear()
                compound.update(**data)
            else:
                compound = Compound(**data)
            return compound
        except Exception as err:
            raise ModiFinderError("Input data is not a valid dictionary. " + str(err)) from err
        

def compound_to_dict(compound):
    """Convert a Compound object to a dictionary"""
    compound_dict = compound.__dict__
    # make sure nothing is passed by reference
    compound_dict = deepcopy(compound_dict)
    return compound_dict


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
            Spectrum = _get_spectrum_class()
            if use_object:
                spectrum = use_object
                spectrum.clear()
                spectrum.update(**data)
            else:
                spectrum = Spectrum(incoming_data=data)
            return spectrum
        except Exception as err:
            raise ModiFinderError("Input data is not a valid dictionary. " + str(err)) from err

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
            raise ModiFinderError("Input data is not a valid Spectrum object. " + str(err)) from err
    
    # USI
    if isinstance(data, str):
        try:
            data = network.get_data(data)
            data = parse_data_to_universal(data)
            data.update(kwargs)
            data = parse_data_to_universal(data)
            Spectrum = _get_spectrum_class()
            if use_object:
                spectrum = use_object
                spectrum.clear()
                spectrum.update(**data)
            else:
                spectrum = Spectrum(**data)
            return spectrum

        except Exception as err:
            raise ModiFinderError("Input data is not a valid USI string. " + str(err)) from err
    
    # Dictionary
    if isinstance(data, dict):
        if needs_parse:
            data = parse_data_to_universal(data)
        try:
            # add kwargs to data
            data.update(kwargs)
            Spectrum = _get_spectrum_class()
            if use_object:
                spectrum = use_object
                spectrum.clear()
                spectrum.update(**data)
            else:
                spectrum = Spectrum(**data)
            return spectrum
        except Exception as err:
            raise ModiFinderError("Input data is not a valid dictionary. " + str(err)) from err
    
    if isinstance(data, list):
        new_data = dict()
        new_data["mz"] = [x[0] for x in data]
        new_data["intensity"] = [x[1] for x in data]
        try:
            Spectrum = _get_spectrum_class()
            if use_object:
                spectrum = use_object
                spectrum.clear()
                spectrum.update(**new_data)
            else:
                spectrum = Spectrum(**new_data)
            return spectrum
        except Exception as err:
            raise ModiFinderError("Input data is not a valid list. " + str(err)) from err
        
    if isinstance(data, _get_spectrum_class()):
        return deepcopy(data)
    
    raise ModiFinderError("Input data is not a valid object.")


    

def spectrum_to_dict(spectrum):
    """Convert a Spectrum object to a dictionary"""
    spectrum_dict = spectrum.__dict__
    # make sure nothing is passed by reference
    spectrum_dict = deepcopy(spectrum_dict)
    return spectrum_dict
