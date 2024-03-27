#!/usr/bin/env python3

import csv
import yaml


representations = {"F": "fundamental", "AS": "antisymmetric", "S": "symmetric"}
channels = {
    "PS": "pseudoscalar",
    "V": "vector",
    "AV": "axialvector",
    "S": "scalar",
    "T": "tensor",
    "AT": "axialtensor",
}


def get_args():
    from argparse import ArgumentParser, FileType
    from sys import stdout

    parser = ArgumentParser()
    parser.add_argument("ensembles", type=FileType("r"))
    parser.add_argument("--output_file", type=FileType("w"), default=stdout)
    return parser.parse_args()


def get_single_line(key, line, content, result):
    result[key] = float(content[line].strip())


def get_file_data(filename):
    result = {
        key: None
        for key in [
            "mass_value",
            "mass_uncertainty",
            "decayconst_value",
            "decayconst_uncertainty",
        ]
    }
    with open(filename, "r") as f:
        content = f.readlines()

    if len(content) >= 6:
        get_single_line("mass_value", 2, content, result)
        get_single_line("mass_uncertainty", 3, content, result)
    if len(content) >= 8:
        get_single_line("decayconst_value", 0, content, result)
        get_single_line("decayconst_uncertainty", 1, content, result)

    return result


def get_single_channel(ensemble, channel):
    filename = (
        "processed_data/Sp{Nc}/beta{beta}/S{NS}T{NT}B{beta}_m{Rep}{m}/"
        "output_{channel}.txt"
    ).format(channel=channel, **ensemble)

    metadata = {key: ensemble[key] for key in ["Nc", "beta", "NS", "NT"]}
    metadata["representation"] = representations[ensemble["Rep"]]
    metadata["valence_mass"] = ensemble["m"]
    metadata["channel"] = channel

    return {**metadata, **get_file_data(filename)}


def collate_values(ensemble, writer):
    for channel, channel_name in channels.items():
        if not ensemble.get(f"Use{channel}"):
            continue

        writer.writerow(get_single_channel(ensemble, channel_name))


def main():
    args = get_args()
    ensembles = yaml.safe_load(args.ensembles)
    writer = csv.DictWriter(
        args.output_file,
        [
            "Nc",
            "beta",
            "NS",
            "NT",
            "representation",
            "channel",
            "valence_mass",
            "mass_value",
            "mass_uncertainty",
            "decayconst_value",
            "decayconst_uncertainty",
        ],
    )
    writer.writeheader()
    for ensemble in ensembles.values():
        collate_values(ensemble, writer)


if __name__ == "__main__":
    main()
