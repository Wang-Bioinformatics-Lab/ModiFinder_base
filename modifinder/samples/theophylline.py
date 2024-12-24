import modifinder as mf
from modifinder.classes.Spectrum import Spectrum
import json

usi = "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00000424785"
accession = "CCMSLIB00000424785"
data = {
    "Compound_Name": "Theophylline",
    "Ion_Source": "DI-ESI",
    "Compound_Source": "Commercial",
    "Instrument": "qTof",
    "PI": "Norberto Lopes",
    "Data_Collector": "Lucas Marques",
    "Adduct": "M+H",
    "Scan": "-1",
    "Precursor_MZ": "181.071",
    "ExactMass": "180.065",
    "Charge": "1",
    "CAS_Number": "58-55-9",
    "Pubmed_ID": " ",
    "Smiles": "O=C1C2=C(N=CN2)N(C)C(N1C)=O",
    "INCHI": "InChI=1S/C7H8N4O2/c1-10-5-4(8-3-9-5)6(12)11(2)7(10)13/h3H,1-2H3,(H,8,9)",
    "INCHI_AUX": " ",
    "Library_Class": "3",
    "SpectrumID": "CCMSLIB00000424785",
    "Ion_Mode": "Positive",
    "create_time": "2014-12-15 03:51:10.0",
    "task_id": "87db0ca37517456084b25b7c32641d73",
    "user_id": "10mauriz",
    "spectrum_id": "CCMSLIB00000424785",
    "source_file": "f.10mauriz/15_br_pos_H_20eV.mgf;",
    "task": "87db0ca37517456084b25b7c32641d73",
    "scan": "1",
    "ms_level": "2",
    "library_membership": "GNPS-LIBRARY",
    "spectrum_status": "1",
    "peaks_json": "[[124.050133,7196.000000],[125.055031,189.000000],[137.081833,615.000000],[142.058365,715.000000],[181.070923,11303.000000],[182.074173,306.000000],[765.732727,376.000000],[1081.180786,123.000000]]",
    "splash": "splash10-0089-0900000000-e75bd6a75670c201ed82",
    "usi": usi,
    "accession": accession,
    "id": accession
}

peaks = json.loads(data["peaks_json"])

spectrum = Spectrum(
    mz= [x[0] for x in peaks],
    intensity= [x[1] for x in peaks],
    precursor_mz= data["Precursor_MZ"],
    precursor_charge= data["Charge"],
    adduct = data["Adduct"],
    ms_level= data["ms_level"],
    instrument= data["Instrument"],
    ms_mass_analyzer= data["Ion_Source"],
    ms_dissociation_method= data["Ion_Mode"],
    spectrum_id= data["SpectrumID"]
)

compound = mf.Compound(structure=data["Smiles"], spectrum=spectrum, name=data["Compound_Name"], id=accession)