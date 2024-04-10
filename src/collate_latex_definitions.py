#!/usr/bin/env python3

from re import search
from sys import stdout


substitutions = {
    "0": "zero",
    "1": "one",
    "2": "two",
    "3": "three",
    "4": "four",
    "5": "five",
    "6": "six",
    "7": "seven",
    "8": "eight",
    "9": "nine",
    ".": "point",
}


def read_captions(filename):
    with open(filename, "r") as f:
        return {
            line.split()[0]: line.split()[1:] for line in f if len(line.split()) > 1
        }


def read_all_captions(filenames):
    captions = {}
    for filename in filenames:
        captions.update(read_captions(filename))
    return captions


def texname(name):
    normalised_name = name
    for symbol, word in substitutions.items():
        normalised_name = normalised_name.replace(symbol, word)
    if invalid_match := search("([^A-Za-z])", normalised_name):
        message = (
            f"{name} could not be converted into something TeX understands "
            f"as it contains the character '{invalid_match.groups()[0]}'"
        )
        raise ValueError(message)
    return f"\\{normalised_name}"


def write_captions(captions, output_file=stdout):
    for name, value in captions.items():
        print(f"\\newcommand {texname(name)} {{{' '.join(value)}}}", file=output_file)


def parse_args():
    from argparse import ArgumentParser, FileType

    parser = ArgumentParser()
    parser.add_argument("input_file", nargs="+")
    parser.add_argument("--output_tex", default=stdout, type=FileType("w"))
    return parser.parse_args()


def main():
    args = parse_args()
    captions = read_all_captions(args.input_file)
    write_captions(captions, args.output_tex)


if __name__ == "__main__":
    main()
