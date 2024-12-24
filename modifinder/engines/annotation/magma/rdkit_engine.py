from rdkit.Chem import *
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, Descriptors
import modifinder.utilities.general_utils as pars
import re
bondtype2string = {v:k for (k, v,) in Chem.rdchem.BondType.names.items()}

class newclass(object):
    """
    Additional class
    """

    pass

def LogP(mol):
    return Chem.Crippen.MolLogP(mol)



def natoms(mol):
    return mol.GetNumAtoms()



def GetExtendedAtomMass(mol, a):
    atom = mol.GetAtomWithIdx(a)
    return pars.mims[atom.GetSymbol()] + pars.Hmass * (atom.GetNumImplicitHs() + atom.GetNumExplicitHs())



def GetAtomSymbol(mol, a):
    return mol.GetAtomWithIdx(a).GetSymbol()



def GetAtomHs(mol, a):
    atom = mol.GetAtomWithIdx(a)
    return atom.GetNumImplicitHs() + atom.GetNumExplicitHs()



def nbonds(mol):
    return mol.GetNumBonds()



def GetBondAtoms(mol, b):
    bond = mol.GetBondWithIdx(b)
    return [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]



def GetNBonds(mol, a):
    return len(mol.GetAtomWithIdx(a).GetBonds())



def GetBondType(mol, b):
    bond = mol.GetBondWithIdx(b)
    return bondtype2string[bond.GetBondType()]



def MolToInchiKey(mol):
    return AllChem.InchiToInchiKey(AllChem.MolToInchi(mol))



def FragmentToSmiles(mol, atomlist):
    emol = Chem.EditableMol(mol)
    for atom in reversed(range(mol.GetNumAtoms())):
        if atom not in atomlist:
            emol.RemoveAtom(atom)

    frag = emol.GetMol()
    for bond in frag.GetBonds():
        # if bond.GetBondType() == Chem.rdchem.BondType.AROMATIC:
        bond.SetBondType(Chem.rdchem.BondType.UNSPECIFIED)
    return Chem.MolToSmiles(frag)



def GetFormulaProps(mol):
    mim = 0.0
    for a in range(mol.GetNumAtoms()):
        mim += GetExtendedAtomMass(mol, a)

    formula_string = Chem.rdMolDescriptors.CalcMolFormula(mol)
    return (mim, formula_string)



def SmilesToMol(smiles, name = None):
    mol = Chem.MolFromSmiles(smiles)
    mol.SetProp('_Name', name)
    AllChem.Compute2DCoords(mol)
    return mol


def GetFormulaMass(formula):
    weight = 0
    pattern = r'([A-Z][a-z]*)(\d*)'
    matches = re.findall(pattern, formula)  # Find all matches in the formula

    # Create a dictionary to store element symbol and count pairs
    for match in matches:
        element = match[0]
        count = match[1]
        if count:
            weight += int(count) * pars.mims[element]
        else:
            weight += pars.mims[element]
    return weight

def get_adduct_mass(adduct):
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
            weight += GetFormulaMass(subformula[1:])
        else:
            weight -= GetFormulaMass(subformula[1:])
    
    if charge[-1] == '+':
        if len(charge) == 1:
            chargeCount = 1
        else:
            chargeCount = int(charge[:-1])
        weight -= pars.elmass * chargeCount
    else:
        if len(charge) == 1:
            chargeCount = 1
        else:
            chargeCount = int(charge[:-1])
        weight += pars.elmass * chargeCount

    return weight