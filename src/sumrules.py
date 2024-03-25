#!/usr/bin/env python3

from collections import defaultdict
import math
from sys import stdout
from uncertainties import ufloat

Ncs = [4, 6, 8]
reps = ["F", "AS", "S"]
observables = ["masses", "decayconsts"]


def s0(vector_decay_const, vector_mass, axialvector_decay_const, axialvector_mass):
    return 4 * math.pi * (
        vector_decay_const_square / vector_mass_square
        - axialvector_decay_const_square / axialvector_mass_square
    )


def s1(vector_decay_const, axialvector_decay_const, pseudoscalar_decay_const):
    return 1 - (
        axialvector_decay_const_square + pseudoscalar_decay_const_square
    ) / vector_decay_const_square


def s2(vector_decay_const, vector_mass, axialvector_decay_const, axialvector_mass):
    return 1 - axialvector_mass_square * axialvector_decay_const_square / (
        vector_mass_square * vector_decay_const_square
    )

s_functions = {"s0": s0, "s1": s1, "s2": s2}


def get_single_N_rep_data(Nc, rep, channel_observable):
    with open(f"processed_data/Sp{Nc}/{rep}/{channel_observable}.dat", "r") as f:
        if len(lines := f.readlines()) == 0 or len(split_line := lines[0].split()) < 2:
            return
        return ufloat(*map(float, split_line))


def get_finite_N_data():
    channel_observables = [
        "pseudoscalar_decayconsts",
        "vector_masses",
        "vector_decayconsts",
        "axialvector_masses",
        "axialvector_decayconsts"
    ]
    return {
        Nc: {
            rep: {
                channel_observable: get_single_file_data(Nc, rep, channel_observable)
                for channel_observable in channel_observables
            }
            for rep in reps
        }
        for Nc in Ncs
    }


def get_single_large_N_data(rep):
    data = {}
    for observable in observables:
        with open(f"processed_data/largeN/{rep}_{observable}.txt", "r") as f:
            data = {**data, **{
                f"{(split_line := line.split())[0]}":
                ufloat(*map(float, split_line[1:3]))
                for line in f
            }}

    return data

def get_large_N_data():
    return {
        rep: get_single_large_N_data(rep)
        for rep in reps
    }

def get_args():
    from argparse import ArgumentParser, FileType
    from sys import stdout

    parser = ArgumentParser()
    parser.add_argument("--output_file", default=stdout, type=FileType("w"))
    return parser.parse_args()


def print_single(label, s_data, s0_slug="", output_file=stdout):
    line_format = r"{label}, ({rep}) & ${s0:.02uSL}$ {s0_slug} & ${s1:.02uSL}$ & ${s2:.02uSL}$ \\"
    print(
        line_format.format(
            label=label,
            s0_slug=s0_slug,
            **s_data
        ),
        file=output_file
    )


def print_bunches(data, slugs=defaultdict(lambda: ""), output_file=stdout):
    for Nc_label, Nc_data in data.items():
        print(r"\hline", file=output_file)
        for rep, s_data in Nc_data.items():
            print_single(
                f"{Nc_label}, ({rep.lower()})",
                s_data,
                s0_slug=slugs[rep],
                output_file=output_file,
            )


def output_table(finite_N, large_N, extras, output_file=stdout):
    print(r"\begin{tabular}{|c|c|c|c|}", file=output_file)
    print(r"\hline", file=output_file)
    print(r"Theory & $s_0$ & $s_1$ & $s_2$ \\", file=output_file)
    print(r"\hline", file=output_file)
    print_bunches(finite_N, output_file=output_file)
    print_bunches(
        large_N,
        slugs={"F": "$N_c$", "AS": "$N_c^2$", "S": "$N_c^2$"},
        output_file=output_file
    )
    print(r"\hline\hline", file=output_file)
    for label, s_data in extras.items():
        print_single(label, s_data)


def compute_single_s(datum):
    return {
        name: function(**datum)
        for name, function in s_functions.items()
    }


def compute_s(data):
    return {
        rep: compute_single_s(datum)
        for rep, datum in data.items()
    }


def main():
    args = get_args()
    finite_N_squared_data = get_finite_N_data()
    large_N_squared_data = get_large_N_data()

    # SU(3) data in MeV, from hep-ph/0501128.pdf
    su3_squared_data = {
        "pseudoscalar_decayconsts": ufloat(92.4, 0.35) ** 2,
        "vector_decayconsts": ufloat(153.4, 7.2) ** 2,
        "vector_masses": ufloat(775.8, 0.5) ** 2,
        "axialvector_decayconsts": ufloat(152.4, 10.4) ** 2,
        "axialvector_masses": ufloat(1230, 40) ** 2,
    }

    finite_N_results = {
        f"$Sp({Nc})$": compute_s(data)
        for Nc, data in finite_N_squared_data.items()
    }
    large_N_results = {r"$Sp(\infty)$": compute_s(large_N_squared_data)}
    su3_results = {
        r"$SU(3)$, $N_{(\mathrm{f})}=2$ ($m_\pi=139.6\textnormal{ MeV}$)":
        compute_single_s(su3_squared_data)
    }
    output_table(
        finite_N_results,
        large_N_results,
        su3_results,
        output_file=args.output_file
    )


if __name__ == "__main__":
    main()
