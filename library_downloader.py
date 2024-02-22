from prettytable import PrettyTable
from urllib.request import urlopen
import matplotlib.pyplot as plt
from rdkit import Chem
from tqdm import tqdm
import pandas as pd
import pickle
from . import utils_n as utils
import json
import os
import sys

def disable_rdkit_logging():
    """
    Disables RDKit whiny logging.
    """
    import rdkit.rdBase as rkrb
    import rdkit.RDLogger as rkl
    logger = rkl.logger()
    logger.setLevel(rkl.ERROR)
    rkrb.DisableLog('rdApp.error')

def get_gnps_library(url, log = True):
    """
    Downloads the GNPS library from the GNPS website.
    """
    if log:
        print("Downloading GNPS library from: " + url)
    response = urlopen(url)
    data_json = json.loads(response.read())
    if log:
        print("Number of compounds in the library: " + str(len(data_json)))
    return data_json

def json_to_dict(data_json, key = 'spectrum_id'):
    """
    Converts the GNPS library JSON to a python dictionary.
    """
    data_dict = dict()
    for i in range(len(data_json)):
        data_dict[data_json[i][key]] = data_json[i]
    return data_dict

def calculate_matches(data_dict, weight_threshold = 1500, difference_threshold_rate = 0.3):
    """
    Calculates the matches between the library.
    It is assumed that the library is passed as a dictionary with unique key identifiers.
    Input:
        data_dict: dictionary with the library
        weight_threshold: maximum weight to consider
        difference_threshold_rate: maximum difference between two weights to consider, if 0 then no threshold is applied
    Output:
        matches: dictionary with the matches, [key] = [list of matches] where key is the number of modification sites
        cachedStructures: dictionary with the cached structures (built from SMILES strings using RDKit)
    """
    matches = {}
    cachedStructures = dict()
    data_ids = list(data_dict.keys())
    for i in tqdm(range(len(data_ids))):   
        
        compound1 = data_ids[i]
        try:
            w1 = float(data_dict[compound1]['Precursor_MZ'])
        except:
            print("Error: while extracting weight for compound " + compound1, file=sys.stderr)
            continue
        if w1 > weight_threshold or data_dict[compound1]['Smiles'] == 'N/A' or data_dict[compound1]['Smiles'] == ' ':
            continue
        
        for j in range(len(data_ids)):
            compound2 = data_ids[j]
            if data_dict[compound1]['Adduct'] != data_dict[compound2]['Adduct']:
                continue
            try:
                w2 = float(data_dict[compound2]['Precursor_MZ'])
            except:
                print("Error: while extracting weight for compound " + compound2, file=sys.stderr)
                continue
            if w2 > weight_threshold or data_dict[compound2]['Smiles'] == 'N/A' or data_dict[compound2]['Smiles'] == ' ':
                continue
            
            if difference_threshold_rate > 0 and abs(w1 - w2) > difference_threshold_rate * min(w1, w2):
                continue

            try:
                if compound1 not in cachedStructures:
                    cachedStructures[compound1] = Chem.MolFromSmiles(data_dict[compound1]['Smiles'])
                if compound2 not in cachedStructures:
                    cachedStructures[compound2] = Chem.MolFromSmiles(data_dict[compound2]['Smiles'])

                m1 = cachedStructures[compound1]
                m2 = cachedStructures[compound2]

                if m1.GetNumAtoms() > m2.GetNumAtoms() and m1.HasSubstructMatch(m2):
                    numModificationSites = len(utils.calculateModificationSites(m1, m2))
                    if numModificationSites not in matches:
                        matches[numModificationSites] = set()
                    matches[numModificationSites].add((compound1, compound2))
            except:
                pass
    
    return matches, cachedStructures

def download(url, output, weight_threshold, difference_threshold_rate, library_name = None, log = True):
    """
    download function.
    Input:
        url: URL to download the JSON library from
        output: output directory
        weight_threshold (mz : double): maximum weight to consider
        difference_threshold_rate (double): maximum difference between two weights to consider
        library_name (string): name of the library
    Output:
        data_dict_filtered: dictionary with the library, filtered to only include the compounds that are in matches
        matches: dictionary with the matches, [key] = [list of matches] where key is the number of modification sites
        cachedStructures: dictionary with the cached structures (built from SMILES strings using RDKit)
    """
    disable_rdkit_logging()

    # if json data is already downloaded, load it
    try:
        with open(os.path.join(output, "jsons", library_name + ".json"), "r") as f:
            data_json = json.load(f)
        if log:
            print("GNPS library already downloaded.")
            print("Number of compounds in the library: " + str(len(data_json)))
    except:
         data_json = get_gnps_library(url, log)
    
    data_dict = json_to_dict(data_json)
    matches, cachedStructures = calculate_matches(data_dict, weight_threshold, difference_threshold_rate)

    # get compound ids that are in matches
    unique_ids = set()
    for key in matches:
        for match in matches[key]:
            unique_ids.add(match[0])
            unique_ids.add(match[1])

    # filter data_dict to only include unique_ids
    data_dict_filtered = dict()
    for key in data_dict:
        if key in unique_ids:
            data_dict_filtered[key] = data_dict[key]
    
    # filter cachedStructures to only include unique_ids
    cachedStructures_filtered = dict()
    for key in cachedStructures:
        if key in unique_ids:
            cachedStructures_filtered[key] = cachedStructures[key]

    print("Number of compounds in the library: " + str(len(data_dict_filtered)))
    t = PrettyTable(['# modification sites', '# matches'])
    for key in matches:
        t.add_row([key, len(matches[key])])
    print(t)
    print("Number of unique ids in matches: " + str(len(unique_ids)))
    
    return data_dict_filtered, matches, cachedStructures_filtered

def calculate_helpers(matches, data_dict_filtered, output, max_num_modifications_allowed = 1):
    """
    Calculates the helpers for each compound.
    Input:
        matches: dataframe of matches
        max_num_modifications_allowed: maximum number of modification sites allowed
    Output:
        helpers: dictionary with the helpers, [key] = [list of helpers] where key is the compound id
    """
    helpers = {}
    for i, row in matches.iterrows():
        t_id_smaller = row["id_smaller"]
        t_id_larger = row["id_bigger"]
        num_modifs = row["num_modif_sites"]
        if t_id_smaller not in helpers:
            helpers[t_id_smaller] = []
        helpers[t_id_smaller].append((t_id_larger, num_modifs))
        if t_id_larger not in helpers:
            helpers[t_id_larger] = []
        helpers[t_id_larger].append((t_id_smaller, num_modifs))
    
    count = []
    # remove duplicates
    for key in helpers:
        helpers[key] = list(set(helpers[key]))
        count.append(len(helpers[key]) - 1)

    # create helpers for each pair (to remove the bigger one)
    for i, row in matches.iterrows():
        if row["num_modif_sites"] > max_num_modifications_allowed:
            continue
        id_smaller = row["id_smaller"]
        id_larger = row["id_bigger"]
        helper_list = list(helpers[id_smaller])
        # remove id_bigger from the list
        to_remove = []
        try:
            for i in range(len(helper_list)):
                # print(i, data_dict_filtered[helper_list[i]]['Precursor_MZ'], data_dict_filtered[id_larger]['Precursor_MZ'])
                if abs(float(data_dict_filtered[helper_list[i][0]]['Precursor_MZ']) - float(data_dict_filtered[id_larger]['Precursor_MZ'])) < 0.5:
                    to_remove.append(helper_list[i][0])
                elif helper_list[i][1] > max_num_modifications_allowed:
                    to_remove.append(helper_list[i][0])
            # print("Number of helpers for {} and {}: {}, and removed: {}".format(id_smaller, id_larger, len(helper_list), len(to_remove)))
        except:
            print("issue with removing", helper_list, to_remove)
            pass
        
        # filter to_remove
        helper_list = [x[0] for x in helper_list if x[0] not in to_remove]
        # if folder does not exist, create it
        if not os.path.exists(os.path.join(output, "helpers")):
            os.makedirs(os.path.join(output, "helpers"))
            # write the list to a jason file
        with open(os.path.join(output, "helpers", "{}_{}.json".format(id_smaller, id_larger)), 'w') as f:
            json.dump(helper_list, f)

def store_cached_values(url, library_name, output, weight_threshold, difference_threshold_rate, log = True):
    """
    Stores the cached values in the output directory.
    Input:
        url: URL to download the JSON library from
        library_name (string): name of the library
        output: output directory
        weight_threshold (mz : double): maximum weight to consider
        difference_threshold_rate (double): maximum difference between two weights to consider
    Output:
        True if successful, False otherwise
    """
    try:
        # create data folder if it does not exist
        if not os.path.exists(output):
            os.makedirs(output)
        
        data_dict_filtered, matches, cachedStructures_filtered = download(url, output, weight_threshold, difference_threshold_rate, library_name, log)


        # print("here1")
        # create matches dataframe
        columns = ['id_bigger', 'id_smaller', 'num_modif_sites', 'weight_bigger', 'weight_smaller', 'difference', "adduct", "weight_threshold", "difference_threshold_rate"]
        matches_df = pd.DataFrame(columns=columns)

        # print("here2")
        for num_modif_sites in matches:
            for match in matches[num_modif_sites]:
                try:
                    weight_bigger = float(data_dict_filtered[match[0]]['Precursor_MZ'])
                    weight_smaller = float(data_dict_filtered[match[1]]['Precursor_MZ'])
                    difference = abs(weight_bigger - weight_smaller)
                    adduct = data_dict_filtered[match[0]]['Adduct']
                    matches_df = pd.concat([matches_df, pd.DataFrame([[match[0], match[1], num_modif_sites, weight_bigger, weight_smaller, difference, adduct, weight_threshold, difference_threshold_rate]], columns=columns)])
                except:
                    pass
        
        # store matches dataframe
        # create directory for matches if it does not exist
        if not os.path.exists(os.path.join(output, "matches")):
            os.makedirs(os.path.join(output, "matches"))
        matches_df.to_csv(os.path.join(output, "matches", library_name + ".csv"), index=False)
        # print("here3")
        def local_store_cached_values(dictionary, name):
            # create directory for dictionary if it does not exist
            if not os.path.exists(os.path.join(output, name)):
                os.makedirs(os.path.join(output, name))
            for key in dictionary:
                with open(os.path.join(output, name, key + ".pkl"), "wb") as f:
                    pickle.dump(dictionary[key], f)
            
        # store compounds_data
        local_store_cached_values(data_dict_filtered, "compounds_data")
        with open(os.path.join(output, "compounds_data", library_name+".pkl"), "wb") as f:
            pickle.dump(data_dict_filtered, f)
        # store cachedStructures
        local_store_cached_values(cachedStructures_filtered, "cached_structures")
        with open(os.path.join(output, "cached_structures", library_name+".pkl"), "wb") as f:
            pickle.dump(cachedStructures_filtered, f)
        
        calculate_helpers(matches_df, data_dict_filtered, output)

        return True
    except:
        print("here")
        return False

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='library downloader')
    parser.add_argument('--url', type=str, default="https://gnps-external.ucsd.edu/gnpslibrary/GNPS-LIBRARY.json", help='URL to download the JSON library from')
    parser.add_argument('--output', type=str, default="Data/GNPS_Library/", help='Output path')
    parser.add_argument('--weight_threshold', type=float, default=500, help='Maximum weight to consider')
    parser.add_argument('--difference_threshold_rate', type=float, default=0.5, help='Maximum difference between two weights to consider')
    args = parser.parse_args()
    download(**vars(args))