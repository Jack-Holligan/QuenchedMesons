#!/usr/bin/env python3

from pathlib import Path
import os
import yaml


def generate_wls(Nc, beta, slug):
    with open("metadata/ensembles.yaml", "r") as f:
        ensembles = yaml.safe_load(f)
    params = ensembles[f"Nc{Nc}_{slug}"]

    input_dir = f"{os.getcwd()}/processed_data/Sp{Nc}/beta{beta}/{slug}"
    output_dir = f"{os.getcwd()}/processed_data/Sp{Nc}/beta{beta}/{slug}"

    with open("src/fit_correlation_function.wls", "r") as f:
        script = f.read()

    for param_name, param_value in params.items():
        script = script.replace(f"_SED_{param_name.upper()}_", str(param_value))
    script = script.replace("_SED_BASEPATH_", input_dir)
    script = script.replace("_SED_OUTPUTPATH_", output_dir)
    return script


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("Nc")
    parser.add_argument("beta")
    parser.add_argument("slug")
    args = parser.parse_args()
    print(generate_wls(args.Nc, args.beta, args.slug))


if __name__ == "__main__":
    main()
