import modifinder as mf
from modifinder.classes.Spectrum import Spectrum
import json

usi = "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00005435812"
accession = "CCMSLIB00005435812"
data = {
    "Compound_Name": "Caffeine",
    "Ion_Source": "LC-ESI",
    "Compound_Source": "Isolated",
    "Instrument": "qTof",
    "PI": "Le Pogam",
    "Data_Collector": "Turpin",
    "Adduct": "M+H",
    "Scan": "-1",
    "Precursor_MZ": "195.087",
    "ExactMass": "194.08",
    "Charge": "1",
    "CAS_Number": "95789-13-2",
    "Pubmed_ID": " ",
    "Smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "INCHI": "1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
    "INCHI_AUX": " ",
    "Library_Class": "3",
    "SpectrumID": "CCMSLIB00005435812",
    "Ion_Mode": "Positive",
    "create_time": "2019-05-03 07:12:03.0",
    "task_id": "7420f79691a94190be1ade8732b33d58",
    "user_id": "pilepoga",
    "spectrum_id": "CCMSLIB00005435812",
    "source_file": "f.pilepogam/Victor_Caffeine/Cafeine2.mgf;",
    "task": "7420f79691a94190be1ade8732b33d58",
    "scan": "1",
    "ms_level": "2",
    "library_membership": "GNPS-LIBRARY",
    "spectrum_status": "1",
    "peaks_json": "[[110.066925,38499.089844],[138.060638,412152.093750],[195.079575,6894530.000000],[195.200180,480874.812500],[196.082092,43027.628906]]",
    "splash": "null-null-null-null",
    "canAdmin": 0,
    "usi": usi,
    "accession": accession,
    "id": accession,
    "spectrum_tags": [
        {
        "tag_database_url": "N/A",
        "tag_task_id": "7f0c0571bc594978a51422ddd5cc25ff",
        "tag_desc": "Food:Other",
        "tag_type": "Sample Type Detected",
        "tag_database": "N/A"
        }
    ]
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