"""
General utility functions
"""
from pyteomics import mgf
import pandas as pd
import numpy as np
import json
import re
import modifinder.utilities.gnps_types as gt
import copy

def is_shifted(val1:float, val2:float, ppm:float=None, mz_tol:float=None) -> bool:
    """
    Determine if two values differ by more than a specified tolerance.
    
    The function checks if the absolute difference between two values exceeds either a given parts per million (ppm) value or a mass/charge (m/z) tolerance. 
    If only one of ppm or mz_tol is provided, the function uses that value for comparison. 
    If both are provided, the function checks both conditions and returns True if either condition is satisfied.
    
    Parameters:
        :val1 (float): The first value to compare.
        :val2 (float): The second value to compare.
        :ppm (float, optional): The parts per million tolerance. Default is None.
        :mz_tol (float, optional): The m/z tolerance. Default is None.
    
    Returns:
        :bool: True if the values differ by more than the specified tolerance, False otherwise.
    
    Raises:
        :ValueError: If neither ppm nor mz_tol is provided.
    """
    diff = abs(val1 - val2)
    if ppm is None:
        if mz_tol is None:
            raise ValueError("Either ppm or mz_tol must be provided")
        return diff > mz_tol
    else:
        if mz_tol is None:
            return diff > max(val1, val2) * ppm / 1e6
        else:
            return diff > mz_tol or diff > max(val1, val2) * ppm / 1e6


def read_mgf(mgf_path: str) -> pd.DataFrame:
    """
    Read an MGF file into a pandas DataFrame
    
    Parameters
    ----------
        mgf_path : str
            The path to the MGF file to be read
    
    Returns
    -------
        pd.DataFrame
                pandas DataFrame with columns as metadata and 'spectrum' as the m/z and intensity values
    """

    msms_df = []
    with mgf.MGF(mgf_path) as reader:
        for spectrum in reader:
            try:
                d = spectrum['params']
                d['spectrum'] = np.array([spectrum['m/z array'],
                                        spectrum['intensity array']])
                if 'precursor_mz' not in d:
                    d['precursor_mz'] = d['pepmass'][0]
                else:
                    d['precursor_mz'] = float(d['precursor_mz'])
                msms_df.append(d)
            except Exception as e:
                print(e)

    msms_df = pd.DataFrame(msms_df)
    if 'precursor_mz' in msms_df.columns and 'scans' in msms_df.columns:
        msms_df['precursor_mz'] = msms_df['precursor_mz'].astype(float)
        msms_df['scans'] = msms_df['scans'].astype(int)
    return msms_df


def write_mgf(msms_df: pd.DataFrame, mgf_path: str):
    """
    Writes a pandas DataFrame to an MGF file

    Parameters
    ----------
        msms_df : pd.DataFrame
            pandas DataFrame with column 'spectrum' and other columns as metadata
        mgf_path : str
            Path to write the MGF file
    """

    specs = []
    for i, row in msms_df.iterrows():
        spectrum = {
            'params': row.drop('spectrum').to_dict(),
            'm/z array': row['spectrum'][0],
            'intensity array': row['spectrum'][1]
        }
        specs.append(spectrum)
    with open(mgf_path, 'w') as out:
        mgf.write(specs, out)

mims = {'H': 1.0078250321,
 'He': 3.016029,
 'Li': 6.015122,
 'Be': 9.012182,
 'B': 10.012937,
 'C': 12.000000,
 'N': 14.0030740052,
 'O': 15.9949146221,
 'F': 18.9984032,
 'Ne': 19.992440,
 'Na': 22.9897692809,
 'Mg': 23.985042,
 'Al': 26.981538,
 'Si': 27.976927,
 'P': 30.97376151,
 'S': 31.97207069,
 'Cl': 34.96885271,
 'Ar': 35.967546,
 'K': 38.96370668,
 'Ca': 39.962591,
 'Sc': 44.955910,
 'Ti': 45.952629,
 'V': 49.947163,
 'Cr': 49.946050,
 'Mn': 54.938050,
 'Fe': 53.939615,
 'Co': 58.933200,
 'Ni': 57.935348,
 'Cu': 62.929601,
 'Zn': 63.929147,
 'Ga': 68.925581,
 'Ge': 69.924250,
 'As': 74.921596,
 'Se': 73.922477,
 'Br': 78.9183376,
 'Kr': 77.920386,
 'Rb': 84.911789,
 'Sr': 83.913425,
 'Y': 88.905848,
 'Zr': 89.904704,
 'Nb': 92.906378,
 'Mo': 91.906810,
 'Tc': 97.907216,
 'Ru': 95.907598,
 'Rh': 102.905504,
 'Pd': 101.905608,
 'Ag': 106.905093,
 'Cd': 105.906458,
 'In': 112.904061,
 'Sn': 111.904821,
 'Sb': 120.903818,
 'Te': 119.904020,
 'I': 126.904468,
 'Xe': 123.905896,
 'Cs': 132.905447,
 'Ba': 129.906310,
 'La': 137.907107,
 'Ce': 135.907144,
 'Pr': 140.907648,
 'Nd': 141.907719,
 'Pm': 144.912744,
 'Sm': 143.911995,
 'Eu': 150.919846,
 'Gd': 151.919788,
 'Tb': 158.925343,
 'Dy': 155.924278,
 'Ho': 164.930319,
 'Er': 161.928775,
 'Tm': 168.934211,
 'Yb': 167.933894,
 'Lu': 174.940768,
 'Hf': 173.940040,
 'Ta': 179.947466,
 'W': 179.946706,
 'Re': 184.952956,
 'Os': 183.952491,
 'Ir': 190.960591,
 'Pt': 189.959930,
 'Au': 196.966552,
 'Hg': 195.965815,
 'Tl': 202.972329,
 'Pb': 203.973029,
 'Bi': 208.980383}

Hmass = mims['H']
elmass = 0.0005486
ionmasses = {1: {'+H': mims['H'],
     '+NH4': mims['N'] + 4 * mims['H'],
     '+Na': mims['Na'],
     '-OH': -(mims['O'] + mims['H']),
     '+K': mims['K']},
 -1: {'-H': -mims['H'],
      '+Cl': mims['Cl']}}


def get_formula_mass(formula):
    """
    Return the mass of a chemical formula.

    Parameters
    ----------
    formula : str
        The chemical formula.
        
    Returns
    -------
    float
        The mass of the chemical formula.
    """
    weight = 0
    pattern = r'([A-Z][a-z]*)(\d*)'
    matches = re.findall(pattern, formula)  # Find all matches in the formula

    # Create a dictionary to store element symbol and count pairs
    for match in matches:
        element = match[0]
        count = match[1]
        if count:
            weight += int(count) * mims[element]
        else:
            weight += mims[element]
    return weight

def get_adduct_mass(adduct):
    """
    Get the mass of an adduct.

    Parameters
    ----------
    adduct : str
        The adduct to get the mass of. 
        
        The supported format is [aM±bA]c± where a and b and c are numbers, M is the molecule, A is the adduct, and c is the charge.
    
    Returns
    -------
    float
        The mass of the adduct.
        
    Raises
    ------
    ValueError
        If the adduct format is not accepted
    """
    weight = 0
    # remove spaces
    adduct = adduct.replace(' ', '')
    acceptedAdductsFormat = re.compile(r'\[M(?:\+[A-Za-z0-9]+|\-[A-Za-z0-9]+)*\][0-9]*[+-]')
    if not acceptedAdductsFormat.match(adduct):
        raise ValueError('Adduct format not accepted')

    charge = adduct.split(']')[1]
    remaining = adduct.split(']')[0]
    remaining = remaining.replace('[','')
    
    regexPattern = re.compile(r'\+[A-Za-z0-9]+|\-[A-Za-z0-9]+')
    subformulas = regexPattern.findall(remaining)
    for subformula in subformulas:
        if subformula[0] == '+':
            weight += get_formula_mass(subformula[1:])
        else:
            weight -= get_formula_mass(subformula[1:])
    
    if charge[-1] == '+':
        if len(charge) == 1:
            chargeCount = 1
        else:
            chargeCount = int(charge[:-1])
        weight -= elmass * chargeCount
    else:
        if len(charge) == 1:
            chargeCount = 1
        else:
            chargeCount = int(charge[:-1])
        weight += elmass * chargeCount

    return weight

def convert_to_universal_key(key: str) -> str:
    """
    Convert different types of keys to universal keys.
    This function standardizes various key names to a universal format. 

    Args:
        :key (str): The key to be converted.
    
    Returns:
        :str: The converted key.
    """
    key = key.lower()
    key = key.replace(" ", "_")
    return gt.gnps_keys_mapping.get(key, key)
    
# TODO: use machine learning prepared data instead of hardcoding
def parse_data_to_universal(data):
    """
    Parse the data to a universal format.

    This function takes a dictionary of data and converts it into a universal format.
    It processes specific keys like "peaks_json" and "Charge" differently, and attempts
    to convert other values to floats. If the conversion to float is successful and the
    key is "Charge", it further converts the value to an integer.

    Args:
        :data (dict): The input data dictionary to be parsed.

    Returns:
        :dict: A dictionary with keys converted to a universal format and values processed
              accordingly.
    """

    res = {}
    for key, value in data.items():
        converted_key = convert_to_universal_key(key)
        if key == "peaks_json":
            res['peaks'] = json.loads(value)
        elif converted_key == "adduct":
            res[converted_key] = gt.adduct_mapping.get(value, value)
        else:
            try:
                if key in ["precursor_charge", "precursor_charge", "ms_level", "scan", "exact_mass"]:
                    value = float(value)
                if key in ["precursor_charge", "charge", "ms_level"]:
                        value = int(value)
            except Exception:
                raise ValueError(f"Could not convert {key} to number")
            res[converted_key] = value
    return res


def entropy(probabilities):
    # if probabilities is not numpy array, convert it to numpy array
    probabilities = np.array(probabilities)
    # if probabilities is not float, convert it to float
    probabilities = probabilities.astype(float)
    
    if len(probabilities) == 0:
        return 1
    if min(probabilities) < 0:
        probabilities = probabilities - min(probabilities)
    regulator = 1e-8
    probabilities = probabilities + regulator
    probabilities = probabilities / np.sum(probabilities)
    H_max = np.log(len(probabilities))
    H = abs(np.sum(probabilities * np.log(probabilities)))
    # print(H, probabilities, H_max, H/H_max)
    return H/H_max


def power_prob(probabilities):
    # copy the probabilities to avoid changing the original
    probabilities2 = copy.deepcopy(probabilities)
    if min(probabilities2) < 0:
        probabilities2 = probabilities2 - min(probabilities2)
    if max(probabilities2) == 0:
        return probabilities2
    # make anythin less than half of the max value zero
    probabilities2[probabilities2 < max(probabilities2) / 2] = 0

    probabilities2 = np.power(probabilities2, 4)
    if sum(probabilities2) == 0:
        return probabilities2
    probabilities2 = probabilities2 / probabilities2.sum()
    return probabilities2


def parse_molecular_formula(formula):
    # Define the regular expression pattern to match element symbols and their counts
    pattern = r'([A-Z][a-z]*)(\d*)'
    matches = re.findall(pattern, formula)  # Find all matches in the formula

    # Create a dictionary to store element symbol and count pairs
    atom_counts = {}
    for match in matches:
        element = match[0]
        element = element.capitalize()
        count = match[1]
        if count:
            count = int(count)
        else:
            count = 1
        atom_counts[element] = count

    return atom_counts


def is_submolecule(sub_formula, target_formula, ignore_H = False):
    # Parse the atom counts of the sub-molecule and target molecule
    sub_atom_counts = parse_molecular_formula(sub_formula)
    target_atom_counts = parse_molecular_formula(target_formula)

    # Check if every atom in sub-molecule is in target molecule and has less or equal count
    for element, count in sub_atom_counts.items():
        if element == 'H' and ignore_H:
            continue
        if element not in target_atom_counts or target_atom_counts[element] < count:
            return False

    return True

def add_adduct_to_formula(formula, adduct):
    if adduct != "M+H":
        raise ValueError("Only H adduct is supported.")
    else:
        adduct = copy.deepcopy(adduct)
        adduct = "H"
    # add one H to the formula
    pattern = r"([A-Z][a-z]*)(\d*)"
    matches = re.findall(pattern, formula)  # Find all matches in the formula
    formula = ""
    Found = False
    for match in matches:
        if match[0] == adduct:
            count = match[1]
            if count:
                count = int(count)
            else:
                count = 1
            count += 1
            formula += adduct + str(count)
            Found = True
        else:
            formula += match[0]
            if match[1]:
                formula += match[1]
    if not Found:
        formula += adduct
    
    return formula

def remove_adduct_from_formula(formula, adduct):
    if adduct != "M+H":
        raise ValueError("Only H adduct is supported.")
    else:
        adduct = copy.deepcopy(adduct)
        adduct = "H"
    
    # remove one H from the formula
    pattern = r'([A-Z][a-z]*)(\d*)'
    matches = re.findall(pattern, formula)  # Find all matches in the formula
    formula = ""
    for match in matches:
        if match[0] == "H":
            count = match[1]
            if count:
                count = int(count)
            else:
                count = 1
            count -= 1
            if count > 1:
                formula += "H" + str(count)
            elif count == 1:
                formula += "H"
        else:
            formula += match[0]
            if match[1]:
                formula += match[1]
    
    return formula

def get_elements(formula):
    elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    elements = [(element, count if count else 1) for element, count in elements]
    return elements

def get_diff(formula1, formula2):
    elements1 = get_elements(formula1)
    elements2 = get_elements(formula2)
    diff = {}
    for element, count in elements1:
        diff[element] = int(count)
    for element, count in elements2:
        if element in diff:
            diff[element] -= int(count)
        else:
            diff[element] = -int(count)
    
    # remove 0s
    diff = {key: value for key, value in diff.items() if value != 0}
    # dict to string
    diff = str(diff)
    return diff

def convert_to_formula(item):
    item = item.replace("{", "")
    item = item.replace("}", "")
    item = item.split(",")
    item = [re.sub(r'[\'\s]', '', i) for i in item]
    item = [i.split(":") for i in item]
    item = [[i[0], int(i[1])] for i in item]
    # item = sorted(item, key=lambda x: x[1], reverse=True)
    item = [i[0] + str(i[1]) for i in item]
    item = "".join(item)
    return item