#!/usr/bin/env python3

import re
import sys


def read_file(filename, row1, row2):
    with open(filename, "r") as f:
        lines = f.readlines()

    if len(lines) < 4:
        return None

    beta, bare_mass = map(float, re.match(".*B([0-9.]+)_m[FAS]+([-0-9.]+)$", filename.split("/")[-2]).groups())
    mass = float(lines[row1])
    mass_error = float(lines[row2])

    return beta, bare_mass, mass, mass_error


def read_mass(filename):
    return read_file(filename, 2, 3)


def read_decayconst(filename):
    return read_file(filename, 0, 1)


def read_w0(filename):
    data = {}
    with open(filename, "r") as f:
        for line in f:
            if len(split_line := line.split()) > 1:
                data[split_line[0]] = float(split_line[1])
    return data


def get_args():
    from argparse import ArgumentParser, FileType
    parser = ArgumentParser()
    parser.add_argument("input_file", nargs="+")
    parser.add_argument("--output_file", type=FileType("w"), default=sys.stdout)
    parser.add_argument("--w0_file", required=True)
    parser.add_argument("--observable", default="mass")
    return parser.parse_args()


def format_output(values, w0):
    output = []
    for beta, bare_mass, value, error in values:
        if beta != w0["beta"]:
            raise ValueError("betas don't match")
        output.append(f"{beta},{bare_mass},{w0['w']},{w0['dw']},{value},{error}")
    return '\n'.join(output)


def main():
    args = get_args()
    readers = {"mass": read_mass, "decayconst": read_decayconst}
    values = [readers[args.observable](input_file) for input_file in args.input_file]
    w0 = read_w0(args.w0_file)
    print(format_output(values, w0), file=args.output_file)


if __name__ == "__main__":
    main()
