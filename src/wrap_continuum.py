#!/usr/bin/env python3

from pathlib import Path
import os
import re
import yaml

all_colours = ["Blue", "Orange", "Darker[Green]", "Red", "Purple"]


def generate_wls(Nc, rep, channel, observable):
    with open("metadata/continuum.yaml", "r") as f:
        metadata = yaml.safe_load(f)
    params = metadata[f"Sp{Nc}"][rep][channel][observable]

    with open("metadata/puregauge.yaml", "r") as f:
        ensembles = yaml.safe_load(f)

    volumes = []
    betas = []
    for ensemble in ensembles:
        sp_tag, slug = ensemble.split("_")
        if sp_tag != f"Sp{Nc}":
            continue
        volumes.append(slug)
        betas.append(float(re.match(".*B(.*)", slug).groups()[0]))

    exclude_volumes = metadata[f"Sp{Nc}"].get("exclude_volumes", [])
    include_volumes = [volume for volume in volumes if volume not in exclude_volumes]

    include_betas = [
        ensembles[f"Sp{Nc}_{volume}"]["beta"] for volume in include_volumes
    ]
    colours = all_colours[: len(include_betas)]

    volumes_str = ", ".join(f'"{volume}"' for volume in include_volumes)
    betas_str = ", ".join(map(str, include_betas))
    colours_str = ", ".join(colours)

    input_dir = f"{os.getcwd()}/processed_data/Sp{Nc}/continuum"
    output_dir = f"{os.getcwd()}/processed_data/Sp{Nc}/continuum"

    with open("src/continuum.wls", "r") as f:
        script = f.read()

    for param_name, param_value in params.items():
        script = script.replace(f"_SED_{param_name.upper()}_", str(param_value))

    script = script.replace("_SED_NC_", str(Nc))
    script = script.replace("_SED_REP_", str(rep))
    script = script.replace("_SED_CHANNEL_", str(channel))
    script = script.replace("_SED_OBSERVABLE_", str(observable))
    script = script.replace("_SED_VOLUMES_", volumes_str)
    script = script.replace("_SED_BETAS_", betas_str)
    script = script.replace("_SED_COLOURS_", colours_str)

    script = script.replace("_SED_BASEPATH_", input_dir)
    script = script.replace("_SED_OUTPUTPATH_", output_dir)
    return script


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("Nc")
    parser.add_argument("rep")
    parser.add_argument("channel")
    parser.add_argument("observable")
    args = parser.parse_args()
    print(generate_wls(args.Nc, args.rep, args.channel, args.observable))


if __name__ == "__main__":
    main()
