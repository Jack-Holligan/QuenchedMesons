#!/usr/bin/env python3

import csv
from functools import partial
from sys import stdout
from uncertainties import ufloat
import yaml


def parse_args():
    from argparse import ArgumentParser, FileType

    parser = ArgumentParser()
    parser.add_argument("ensembles_file")
    parser.add_argument("--output_table", type=FileType("w"), default=stdout)
    parser.add_argument("--output_csv", type=FileType("w"), default=None)
    return parser.parse_args()


def get_w0_from_file(filename):
    with open(filename, "r") as f:
        content = {
            split_line[0]: float(split_line[1])
            for line in f
            if len(split_line := line.split()) > 1
        }
    if "w" not in content or "dw" not in content:
        raise ValueError(f"File {filename} has missing content.")

    return ufloat(content["w"], content["dw"])


def get_w0(ensemble):
    filename = "processed_data/Sp{Nc}/beta{beta}/wflow.dat"
    return {
        **ensemble,
        "w0": get_w0_from_file(filename.format(**ensemble)),
    }


def get_w0s(ensembles):
    return [get_w0(ensemble) for ensemble in ensembles.values()]


def group_by_Nc(w0s):
    return {
        Nc: sorted([w0 for w0 in w0s if w0["Nc"] == Nc], key=lambda w: w["beta"])
        for Nc in sorted(set(w0["Nc"] for w0 in w0s))
    }


def tabulate(w0s, f=stdout):
    print(r"\begin{tabular}{|c|c|c|c|}", file=f)
    print(r"\hline", file=f)
    print(r"$N$ & $\beta$ & $N_s^3 \times N_t$ & $w_0 / a$ \\", file=f)
    print(r"\hline", file=f)
    for Nc, Nc_w0s in group_by_Nc(w0s).items():
        num_rows = len(Nc_w0s)
        row_template = ["", "{beta}", r"${Ns}^3 \times {Nt}$", "${w0:.02uSL}$"]
        contents = [list(map(lambda e: e.format(**w0), row_template)) for w0 in Nc_w0s]
        contents[0][0] = f"\\multirow{{{num_rows}}}*{{{Nc // 2}}}"

        print(r"\hline", file=f)
        print(
            " \\\\\n\\cline{2-4}".join(
                [" & ".join(row_content) for row_content in contents]
            ),
            file=f,
        )
        print("\\\\\n\\hline", file=f)
    print(r"\end{tabular}", file=f)


def write_csv(w0s, output_file):
    writer = csv.DictWriter(
        output_file,
        ["Nc", "beta", "NS", "NT", "w0_value", "w0_uncertainty"],
    )
    writer.writeheader()
    for w0 in w0s:
        writer.writerow(
            {
                "Nc": w0["Nc"],
                "beta": w0["beta"],
                "NS": w0["Ns"],
                "NT": w0["Nt"],
                "w0_value": w0["w0"].nominal_value,
                "w0_uncertainty": w0["w0"].std_dev,
            }
        )


def main():
    args = parse_args()
    with open(args.ensembles_file, "r") as f:
        ensembles = yaml.safe_load(f)
    w0s = get_w0s(ensembles)
    tabulate(w0s, f=args.output_table)
    if args.output_csv:
        write_csv(w0s, args.output_csv)


if __name__ == "__main__":
    main()
