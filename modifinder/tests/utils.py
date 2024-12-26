import unittest


def compare_compounds(compound1, compound2):
    """Compares two compounds and returns True if they are the same."""
    for key in ["adduct_mass", "is_known", "name"]:
        if getattr(compound1, key) != getattr(compound2, key):
            print(f"Key {key} does not match", getattr(compound1, key), getattr(compound2, key))
            return False
    
    if compound1.structure is not None and compound2.structure is not None:
        if not compound1.structure.HasSubstructMatch(compound2.structure):
            print("Structures do not match")
            return False
        if not compound2.structure.HasSubstructMatch(compound1.structure):
            print("Structures do not match")
            return False
    elif compound1.structure is not None or compound2.structure is not None:
        print("One compound has a structure and the other does not")
        return False

    if compound1.spectrum is not None and compound2.spectrum is not None:
        if compound1.spectrum.mz != compound2.spectrum.mz:
            print("MZ does not match", compound1.spectrum.mz, compound2.spectrum.mz)
            return False
        if compound1.spectrum.intensity != compound2.spectrum.intensity:
            print("Intensity does not match")
            return False
        if compound1.spectrum.precursor_mz != compound2.spectrum.precursor_mz:
            print("Precursor MZ does not match")
            return False
        if compound1.spectrum.precursor_charge != compound2.spectrum.precursor_charge:
            print("Precursor charge does not match")
            return False
        if compound1.spectrum.adduct != compound2.spectrum.adduct:
            print("Adduct does not match")
            return False
    elif compound1.spectrum is not None or compound2.spectrum is not None:
        print("One compound has a spectrum and the other does not")
        return False
    
    return True