import os
from pathlib import Path
import subprocess
import pandas as pd
from diffexpr.py_deseq import py_DESeq2
# Takes already pre-filtered dataset, and sets up input for DESeq2
# Assume looked through # of samples etc and specified correct design



def run_command(args):
    """Run command, transfer stdout/stderr"""
    result = subprocess.run(args)
    try:
        result.check_returncode()
    except subprocess.CalledProcessError as e:
        raise e


def run_deseq2(fitness_dir, experiment_name, sdf, edf, design, r_path, feat_id):
    sdf_path = Path(fitness_dir) / f"{experiment_name}_sdf.csv"
    edf_path = Path(fitness_dir) / f"{experiment_name}_edf.csv"
    sdf.to_csv(sdf_path)
    edf.set_index(feat_id).to_csv(edf_path)
    cmd = f'Rscript {r_path} {sdf_path} {edf_path} {experiment_name} {design} {fitness_dir}'
    print(cmd)
    r = run_command(cmd.split())
    os.remove(sdf_path)
    os.remove(edf_path)


def get_deseq2_results(outDir, experiment_name):
    fitness_files = [f for f in outDir.glob(f"{experiment_name}*results*csv")]
    fitness_df = (pd.concat([pd.read_csv(f, sep=' ').assign(day=f.stem.split("-")[-1]) for f in fitness_files])
                    .assign(experiment=experiment_name)
                    .reset_index())
    return fitness_df

