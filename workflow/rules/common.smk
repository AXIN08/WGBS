
import pandas as pd
import os
import argparse
import sys

configfile: "/home/liuxin/data_pasteur/12_epigenome/WGBS/workflow/config/config.yaml"

def parse_samples(samples_tsv):
    return pd.read_csv(samples_tsv, sep='\t').set_index("id", drop=False)

def get_path(sample_df, wildcards, col):
    return sample_df.loc[[wildcards.sample], col].dropna().values[0]

samples_df = parse_samples(config["sample"])


OUTPUT_DIR = config["outdir"]
REFERENCE_DIR = config["ref"]
SAMPLES = samples_df['id'].tolist()

