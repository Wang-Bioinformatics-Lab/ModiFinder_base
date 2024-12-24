import argparse
import depricated.ModificationSiteLocator as SiteLocator
from depricated.Compound import Compound
import modifinder.visualizer as visualizer
import modifinder.utils as utils

import os
import json
import pickle
import requests
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import Draw
from copy import deepcopy
import matplotlib.pyplot as plt
from IPython.display import SVG

def sitelocalize(usi1, usi2, smiles1, smiles2 = None, show_img=True):

    # make the compounds from the usi and smiles
    mainCompound = Compound(usi1, smiles1)
    modCompound = Compound(usi2, smiles2)

    modsite = SiteLocator.ModificationSiteLocator(mainCompound, modCompound)
    likelihoods = modsite.generate_probabilities()
    svg = visualizer.highlightScores(mainCompound.structure, likelihoods)
    
    if show_img:
        SVG(svg)
    
    result = {}
    result["# matched peaks"] = len(modsite.matched_peaks)
    result['# unshifted peaks'] = len(modsite.unshifted)
    result['# shifted peaks'] = len(modsite.shifted)
    result["usi1"] = usi1
    result["usi2"] = usi2
    result["structure1"] = smiles1
    
    result['image'] = svg

    return result

def main(s1, s2, smiles, output_file, output_image_folder, show_img=True):
    result = sitelocalize(s1, s2, smiles, show_img=show_img)

    if output_file.endswith('.xlsx'):
        visualizer.table_to_xlsx(result, output_file)
    elif output_file.endswith('.tsv'):
        # writing out the tsv file
        visualizer.table_to_tsv(result, output_file, output_image_folder)


if __name__ == '__main__':
    # arguments
    parser = argparse.ArgumentParser(description='Substructure Assignment')
    parser.add_argument('--s1', type=str, default='mzspec:GNPS:TASK-e31d8bbfb26c4c3e8f159f06e37d65bb-spectra/specs_ms.mgf:scan:3189', help='gnps accepted usi for known molecule')
    parser.add_argument('--s2', type=str, default='mzspec:GNPS:TASK-e31d8bbfb26c4c3e8f159f06e37d65bb-spectra/specs_ms.mgf:scan:3399', help='gnps accepted usi for modified molecule')
    parser.add_argument('--smiles', type=str, default='OC(c1c(CCCCC)cc(OC(c2c(OC)cc(OC)cc2CCCCC)=O)cc1OC)=O', help='smiles for the known molecule')
    parser.add_argument('--output_file', type=str, default="output.xlsx", help='output file')
    parser.add_argument('--output_image_folder', type=str, default=".", help='output img folder')
    parser.add_argument('--show_img', type=bool, default=True, help='show final image')
    args = parser.parse_args()

    main(args.s1, args.s2, args.smiles, args.output_file, args.output_image_folder)