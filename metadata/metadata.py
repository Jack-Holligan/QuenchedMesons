#!/usr/bin/env python3

from collections import defaultdict
import os
import yaml


channel_names = {
    "PS": "pseudoscalar",
    "S": "scalar",
    "V": "vector",
    "AV": "axialvector",
    "T": "tensor",
    "AT": "axialtensor",
}


def freeze_dd(dd):
    for key, element in dd.items():
        if isinstance(element, defaultdict):
            freeze_dd(element)
        elif isinstance(element, set):
            dd[key] = sorted(element)

    dd.default_factory = None


with open(os.environ.get("SP2N_METADATA_FILE", "metadata/ensembles.yaml"), "r") as f:
    metadata = yaml.safe_load(f)

betas = defaultdict(set)
bare_masses = defaultdict(
    lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
)
ensembles = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))

for label, ens in metadata.items():
    if not ens:
        continue

    betas[ens["Nc"]].add(ens["beta"])

    shortslug = lambda e: f"S{e['NS']}T{e['NT']}B{e['beta']}"
    for channel in "AT", "AV", "PS", "S", "T", "V":
        if ens[f"Use{channel}"]:
            bare_masses[ens["Nc"]][shortslug(ens)][ens["Rep"]][
                channel_names[channel]
            ].add(ens["m"])
            ensembles[ens["Nc"]][ens["Rep"]][channel_names[channel]].add(shortslug(ens))

freeze_dd(betas)
freeze_dd(bare_masses)
freeze_dd(ensembles)
