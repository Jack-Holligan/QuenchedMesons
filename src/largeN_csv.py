#!/usr/bin/env python3

import csv
import itertools


def get_args():
    from argparse import ArgumentParser, FileType

    parser = ArgumentParser()
    parser.add_argument("--output_csv", required=True, type=FileType("w"))
    return parser.parse_args()


def to_number_or_name(element):
    channels = {
        "PS": "pseudoscalar",
        "V": "vector",
        "AV": "axialvector",
        "S": "scalar",
        "T": "tensor",
        "AT": "axialtensor",
    }
    if element == "--":
        return None
    elif element in channels:
        return channels[element]
    else:
        return float(element)


def get_file_data(representation, observable):
    plurals = {
        "mass_hat_squared": "masses",
        "decayconst_hat_squared_over_Nc": "decayconsts",
    }
    representations = {"AS": "antisymmetric", "F": "fundamental", "S": "symmetric"}
    filename = f"processed_data/largeN/{representation}_{plurals[observable]}.txt"
    try:
        with open(filename, "r") as f:
            lines = f.readlines()
    except FileNotFoundError:
        return

    data = []
    for line in lines:
        datum = {
            "representation": representations[representation],
            "observable": observable,
        }
        for label, element in zip(
            [
                "channel",
                "large_N_value",
                "large_N_uncertainty",
                "chisquare",
                "Delta_value",
                "Delta_uncertainty",
            ],
            line.split(),
        ):
            datum[label] = to_number_or_name(element)

        data.append(datum)

    return data


def get_all_data():
    representations = "F", "AS", "S"
    observables = "mass_hat_squared", "decayconst_hat_squared_over_Nc"
    return [
        datum
        for representation, observable in itertools.product(
            representations, observables
        )
        for datum in get_file_data(representation, observable)
    ]


def write_csv(data, output_file):
    # Large-N values and Deltas as defined in Eqs. (19)--(22) of 2312.08465
    writer = csv.DictWriter(
        output_file,
        [
            "representation",
            "observable",
            "channel",
            "large_N_value",
            "large_N_uncertainty",
            "Delta_value",
            "Delta_uncertainty",
            "chisquare",
        ],
    )
    writer.writeheader()
    for datum in data:
        if datum:
            writer.writerow(datum)


def main():
    args = get_args()
    data = get_all_data()
    write_csv(data, args.output_csv)


if __name__ == "__main__":
    main()
