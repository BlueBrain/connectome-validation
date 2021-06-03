#!/usr/bin/env python

import sys
import io
import os
import pandas
import subprocess
import shutil
import yaml

import bluepy.v2 as bluepy


extra_args = ["-n", "250", "--assume-syns-bouton", "1.2"]


def reformat(results):
    idx = pandas.Index(results['mtype'])
    series_mn = pandas.Series(results['mean'].values, index=idx)
    series_sd = pandas.Series(results['std'].values, index=idx)
    return series_mn, series_sd


def get_boutons_mdl(circ_path, target_name):
    target_args = ["--sample-target", target_name]
    if shutil.which("connectome-stats") is None:
        print("""This script uses the 'connectome-stats' script, which is part of 'connectome-tools'.
        Please load the connectome-tools module!""")
    print("Calling connectome-stats. This might take a while...")
    output = subprocess.check_output(["connectome-stats", "bouton-density", "--short"] +
                                     extra_args + target_args + [circ_path])
    results = pandas.read_csv(io.StringIO(str(output, "utf-8")), sep='\s+')

    return reformat(results)


def get_bioname(circ_path):
    with open(circ_path, "r") as fid:
        for ln in fid.readlines():
            splt = ln.strip().split()
            if len(splt) == 2:
                if splt[0] == "BioName":
                    return splt[1]
    return os.path.join(os.path.split(circ_path)[0], 'bioname')


def get_bio_data_path(circ_path):
    fn = os.path.join(get_bioname(circ_path), "s2f.yaml")
    with open(fn, "r") as fid:
        data = yaml.load(fid, Loader=yaml.SafeLoader)
    return data["estimate_individual_bouton_reduction"]["bio_data"]


def get_bio_data(circ_path, bio_data=None):
    if bio_data is None:
        try:
            bio_data = get_bio_data_path(circ_path)
            print("""Will use {0} as biological reference dataset. If that is wrong, provide it explicitly:
            {1} CircuitConfig /path/to/reference/dataset.tsv
            Formatting of reference data follows s2f-recipe
            (see documentation of connectome-tools).""".format(bio_data, __file__))
        except:
            print("""Lookup of biological reference dataset failed. Please provide it explicitly using:
            {0} CircuitConfig /path/to/reference/dataset.tsv
            Formatting of reference data follows s2f-recipe
            (see documentation of connectome-tools).""".format(__file__))
    data = pandas.read_csv(bio_data, sep='\s+')
    data.columns = [x.strip() for x in data.columns]
    return reformat(data)


def validate(circ_path, target_name, out_fn, bio_data=None):
    mn_bio, sd_bio = get_bio_data(circ_path, bio_data=bio_data)
    mn_bio.name = "Mean (ref.)"
    sd_bio.name = "Std (ref.)"
    mn_mdl, sd_mdl = get_boutons_mdl(circ_path, target_name)
    mn_mdl.name = "Mean (data)"
    sd_mdl.name = "Std (data)"
    err = 100 * (mn_mdl - mn_bio) / (0.5 * (sd_mdl + sd_bio))
    err.name = "Error"
    result = pandas.concat([mn_mdl, sd_mdl, mn_bio, sd_bio, err], axis=1)
    result.to_csv(out_fn)


def main():
    if len(sys.argv) < 4:
        print("""Usage: {0} CircuitConfig target OutputFilename.csv (path/to/biological/reference/data.tsv)
        If no biological reference is provided, will try to look it up from the circuit bioname. 
        But that might fail, if circuit files not at the standardized location.
        If provided, the biological reference must be in the same format used by s2f-recipe
        (see connectome-tools documentation).""".format(__file__))
        sys.exit(2)
    if len(sys.argv) == 4:
        print("""No biological reference dataset provided!
        Will try to look up relevant data from the circuit bioname""")
    validate(*sys.argv[1:])


if __name__ == "__main__":
    main()

