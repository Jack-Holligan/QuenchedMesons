#!/usr/bin/env python3

from collections import defaultdict
import sys
import yaml
from random import randint

get_seed = False

metadata = defaultdict(
    lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
)


def stripreader(filename):
    for line in open(filename, "r"):
        yield line.strip()


def freeze_dd(dd):
    result = {}
    for key, element in dd.items():
        if isinstance(element, defaultdict):
            result[key] = freeze_dd(element)
        else:
            result[key] = element
    return result


name_map = {
    "Pseudoscalar": "pseudoscalar",
    "Vector": "vector",
    "Axial-vector": "axialvector",
    "Axial-tensor": "axialtensor",
    "Scalar": "scalar",
    "Tensor": "tensor",
    "mass": "masses",
    "decay": "decayconsts",
    "AS": "AS",
    "S": "S",
    "F": "F",
}


def get_metadata_node(Nc, rep, channel, observable):
    return metadata[f"Sp{Nc}"][name_map[rep]][name_map[channel]][name_map[observable]]


for line in stripreader("notebooks_cat"):
    if (
        line.endswith('"Section",')
        and len(split_line := line.split('"')[1].split()) > 1
    ):
        channel, observable = split_line
    elif line.startswith('RowBox[{"Nc", ":=",'):
        Nc = int(line.split('"')[5])
    elif line == 'RowBox[{"SeedRandom", "[",':
        get_seed = True
    elif line.startswith('RowBox[{"Rep", ":="'):
        rep = line.split('"')[6].strip(r"\<>")
    elif line.startswith('RowBox[{"minmassPS", "=",'):
        minmassPS = float(line.split('"')[5])
        get_metadata_node(Nc, rep, channel, observable)["minmassPS"] = minmassPS
    elif line.startswith('RowBox[{"maxmassPS", "=",'):
        minmassPS = float(line.split('"')[5])
        get_metadata_node(Nc, rep, channel, observable)["maxmassPS"] = minmassPS
    elif get_seed:
        get_seed = False
        seed = int(line.split('"')[1])
        if Nc == 8:
            seed = randint(0, 32768)
        get_metadata_node(Nc, rep, channel, observable)["seed"] = seed

metadata = freeze_dd(metadata)
with open("continuum.yaml", "w") as f:
    print(yaml.dump(metadata), file=f)
