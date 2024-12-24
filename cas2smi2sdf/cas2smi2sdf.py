import numpy as np
import pandas as pd

from pubchempy import get_properties
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.SaltRemover import SaltRemover
from openbabel import pybel

import os
from pathlib import Path
from argparse import ArgumentParser
import requests
from datetime import datetime
import pprint

"""
args:
    mode: cas2smi or smi2sdf or cas2sdf
    --url: slack notification API (if you want to use)

    if cas2smi:
        --cas_csv: path of the csv with CAS (!CONFIRM CSV NAME IS YOUR PROJECT NAME!)
        --cas_col: column name of the CAS

    if smi2sdf:
        --smi_csv: path of the csv with SMILES
        --smi_col: column name of the SMILES
        --sdf_dir: directory path of the sdfs (don't have to create dir at first)

    if cas2sdf:
        --cas_csv: path of the csv with CAS (!CONFIRM CSV NAME IS YOUR PROJECT NAME!)
        --cas_col: column name of the CAS
        --sdf_dir: directory path of the sdfs (don't have to create dir at first)
"""


# cas to smiles
def cas2smi(args):
    """_summary_

    Parameters
    ----------
    args

    Returns
    -------
    df_smiles : pd.DataFrame
        dataframe contains SMILES and other information
    """

    slack = True if args.url != "None" else False

    cas = pd.read_csv(args.cas_csv)
    cas = cas[args.cas_col]
    print(cas)
    cas.dropna(inplace=True)
    cas.drop_duplicates(inplace=True)
    print(cas)
    length = len(cas)
    df_smiles = pd.DataFrame([], columns=["CID", "CAS", "IsomericSMILES", "XLogP", "MolecularWeight", "MolecularFormula", "IUPACName"])
    properties = ["IsomericSMILES", "XLogP", "MolecularWeight", "MolecularFormula", "IUPACName"]
    
    if slack:
        payload = {"text": f"cas2smi start: {length}"}
        r = requests.post(url=args.url, json=payload)

    # get_properties
    for i, c in enumerate(cas):
        try:
            info = get_properties(properties, c, "name", as_dataframe=True)
        except:
            if slack:
                if (i+1) % 10000 == 0:
                    payload = {"text": f"cas2smi {i+1}/{length} done."}
                    r = requests.post(url=args.url, json=payload)
            continue
        
        info["CAS"] = c
        info["CID"] = info.index

        df_smiles = pd.concat([df_smiles, info], axis=0)
        
        if slack:
            if (i+1) % 10000 == 0:
                payload = {"text": f"cas2smi {i+1}/{length} done."}
                r = requests.post(url=args.url, json=payload)
    
    df_smiles.reset_index(inplace=True, drop=True)
    df_smiles.to_csv(os.path.join(os.path.dirname(args.cas_csv), f"{datetime.now().strftime("%Y%m%d")}_{args.proj_name}_SMILES.csv"))
    
    if slack:
        payload = {"text": "cas2smi finished."}
        r = requests.post(url=args.url, json=payload)

    return df_smiles


# split
def salt_remove(smiles):
    """_summary_

    Parameters
    ----------
    smiles : str
        SMILES

    Returns
    -------
    split : str
        SMILES without salts or the largest part of the SMILES
    """
    remover = SaltRemover()
    try:
        s2 = Chem.MolToSmiles(remover.StripMol(Chem.MolFromSmiles(smiles), dontRemoveEverything=True), isomericSmiles=True)
        if "." in s2:
            try:
                mol_frags = Chem.GetMolFrags(Chem.MolFromSmiles(s2), asMols=True)
                largest = None
                largest_size = 0
                for mol in mol_frags:
                    size = mol.GetNumAtoms()
                    if size > largest_size:
                        largest = mol
                        largest_size = size
                s2 = Chem.MolToSmiles(largest)
                return s2
            except:
                s2 = s2.split(".")
                max_index = max(range(len(s2)), key=lambda i: len(s2[i]))
                return s2[max_index]
    except:
        smiles = smiles.split(".")
        max_index = max(range(len(smiles)), key=lambda i: len(smiles[i]))
        return smiles[max_index]
    return s2


# smiles to sdf
def smi2sdf(df_smiles, args):
    """_summary_

    Parameters
    ----------
    df_smiles : pd.DataFrame
        dataframe contains SMILES
    args

    Returns
    -------
    failed : dict
        SMILES which failed to get sdf
    """

    slack = True if args.url != "None" else False

    smiles = df_smiles[args.smi_col]
    
    smiles.apply(salt_remove)
    length = len(smiles)

    if slack:
        payload = {"text": f"smi2sdf start: {length}"}
        r = requests.post(url=args.url, json=payload)

    failed = {}

    for i, smi in enumerate(smiles):
        try:
            ligs = pybel.readstring("smi", smi)
            ligs.localopt()
        except:
            failed[i] = smi
            if slack:
                if (i+1) % 10000 == 0:
                    payload = {"text": f"smi2sdf {i+1}/{length} done."}
                    r = requests.post(url=args.url, json=payload)
            continue

        fileout = os.path.join(args.sdf_dir, f"{i}.sdf")
        ligs.write("sdf", fileout)

        if slack:
            if (i+1) % 10000 == 0:
                payload = {"text": f"smi2sdf {i+1}/{length} done."}
                r = requests.post(url=args.url, json=payload)
    
    if slack:
        payload = {"text": "smi2sdf finished."}
        r = requests.post(url=args.url, json=payload)

    return failed


def main(args):

    if args.mode == "cas2sdf":
        df_smiles = cas2smi(args)
        args.smi_col = "IsomericSMILES"
        os.makedirs(args.sdf_dir, exist_ok=True)
        failed = smi2sdf(df_smiles, args)
        if len(failed) != 0:
            if len(failed) <= 10:
                pprint.pprint(failed, width=1)
            else:
                print("failed smi were more than 10")
            with open(os.path.join(os.path.dirname(args.cas_csv), f"{datetime.now().strftime("%Y%m%d")}_{args.proj_name}_failed.txt"), "w") as f:
                for idx, smi in failed.items():
                    f.write(str(idx) + "  " + str(smi) + "\n")
        else:
            print("no failed smi")

    elif args.mode == "cas2smi":
        _ = cas2smi(args)

    elif args.mode == "smi2sdf":
        os.makedirs(args.sdf_dir, exist_ok=True)
        df_smiles = pd.read_csv(args.smi_csv)
        failed = smi2sdf(df_smiles, args)
        if len(failed) != 0:
            if len(failed) <= 10:
                pprint.pprint(failed, width=1)
            else:
                print("failed smi were more than 10")
            with open(os.path.join(os.path.dirname(args.cas_csv), f"{datetime.now().strftime("%Y%m%d")}_{args.proj_name}_failed.txt"), "w") as f:
                for idx, smi in failed.items():
                    f.write(str(idx) + "  " + str(smi) + "\n")
        else:
            print("no failed smi")
    
    else:
        raise NameError('invalid mode. --mode should be "cas2sdf" or "cas2smi" or "smi2sdf"')

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--mode")
    parser.add_argument("--url")
    parser.add_argument("--cas_csv")
    parser.add_argument("--cas_col")
    parser.add_argument("--smi_csv")
    parser.add_argument("--smi_col")
    parser.add_argument("--sdf_dir")
    args = parser.parse_args()
    if args.cas_csv is not None:
        args.proj_name = Path(args.cas_csv).stem
    else:
        args.proj_name = Path(args.smi_csv).stem
    
    main(args)
