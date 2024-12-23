#!/bin/bash

# activate your python environment
source /opt/pip-env/bin/activate

<< COMMENT
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
COMMENT


python3 $(dirname $0)/cas2smi2sdf.py \
    mode cas2sdf \
    --url None \
    --cas_csv None \
    --cas_col None \
    --smi_csv: None \
    --smi_col None \
    --sdf_dir None