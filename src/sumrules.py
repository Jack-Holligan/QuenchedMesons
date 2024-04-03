#!/usr/bin/env python3

from collections import defaultdict
import csv
import math
from sys import stdout
from uncertainties import ufloat

Ncs = [4, 6, 8]
reps = {"F": "fundamental", "AS": "antisymmetric", "S": "symmetric"}
observables = ["masses", "decayconsts"]
channel_names = {
    "PS": "pseudoscalar",
    "V": "vector",
    "AV": "axialvector",
    "T": "tensor",
    "AT": "axialtensor",
    "S": "scalar",
}


def value_or_none(func):
    def wrapped_function(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except TypeError:
            return None

    return wrapped_function


@value_or_none
def s0(
    vector_decayconsts, vector_masses, axialvector_decayconsts, axialvector_masses, **_
):
    return (
        4
        * math.pi
        * (
            vector_decayconsts / vector_masses
            - axialvector_decayconsts / axialvector_masses
        )
    )


@value_or_none
def s1(vector_decayconsts, axialvector_decayconsts, pseudoscalar_decayconsts, **_):
    return 1 - (axialvector_decayconsts + pseudoscalar_decayconsts) / vector_decayconsts


@value_or_none
def s2(
    vector_decayconsts, vector_masses, axialvector_decayconsts, axialvector_masses, **_
):
    return 1 - axialvector_masses * axialvector_decayconsts / (
        vector_masses * vector_decayconsts
    )


s_functions = {"s0": s0, "s1": s1, "s2": s2}


def get_single_N_rep_data(Nc, rep, channel_observable):
    try:
        with open(
            f"processed_data/Sp{Nc}/continuum/{rep}/{channel_observable}_{rep}_Sp{Nc}.dat",
            "r",
        ) as f:
            if (
                len(lines := f.readlines()) == 0
                or len(split_line := lines[0].split()) < 2
                or "--" in split_line
            ):
                return
    except FileNotFoundError:
        return

    return ufloat(*map(float, split_line))


def get_finite_N_data():
    channel_observables = [
        "pseudoscalar_decayconsts",
        "vector_masses",
        "vector_decayconsts",
        "axialvector_masses",
        "axialvector_decayconsts",
    ]
    return {
        Nc: {
            rep: {
                channel_observable: get_single_N_rep_data(Nc, rep, channel_observable)
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
            for line in f:
                split_line = line.split()
                data[f"{channel_names[split_line[0]]}_{observable}"] = ufloat(
                    *map(float, split_line[1:3])
                )

    return data


def get_large_N_data():
    return {rep: get_single_large_N_data(rep) for rep in reps}


def get_args():
    from argparse import ArgumentParser, FileType
    from sys import stdout

    parser = ArgumentParser()
    parser.add_argument("--output_table", default=stdout, type=FileType("w"))
    parser.add_argument("--output_csv", default=None, type=FileType("w"))
    return parser.parse_args()


def print_single(label, s_data, s0_slug="", output_file=stdout):
    line_format = r"{label}"
    for s_label, s_datum in s_data.items():
        if s_datum is not None:
            line_format += f" & ${{{s_label}:.02uSL}}$"
            if s_label == "s0":
                line_format += f" {{{s_label}_slug}}"
        else:
            line_format += r" & $\cdots$"
    line_format += r" \\"

    print(line_format.format(label=label, s0_slug=s0_slug, **s_data), file=output_file)


def print_bunches(data, slugs=defaultdict(lambda: ""), output_file=stdout):
    for Nc, Nc_data in data.items():
        print(r"\hline", file=output_file)
        for rep, s_data in Nc_data.items():
            print_single(
                f"$Sp({Nc})$, ({reps[rep].lower()})",
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
        output_file=output_file,
    )
    print(r"\hline\hline", file=output_file)
    for label, s_data in extras.items():
        print_single(label, s_data, output_file=output_file)
    print(r"\hline", file=output_file)
    print(r"\end{tabular}", file=output_file)


def compute_single_s(datum):
    return {name: function(**datum) for name, function in s_functions.items()}


def compute_s(data):
    return {
        rep: compute_single_s(datum) if datum else None for rep, datum in data.items()
    }


def split_ufloat(prefix, data):
    return {
        f"{prefix}_value": data.nominal_value if data else None,
        f"{prefix}_uncertainty": data.std_dev if data else None,
    }


def write_csv_row(group_family, Nc, representation, Nf, data, writer):
    to_write = {
        "group_family": group_family,
        "Nc": Nc,
        "representation": reps[representation],
        "Nf": Nf,
    }
    for s in "s0", "s1", "s2":
        to_write.update(**split_ufloat(s, data[s]))

    writer.writerow(to_write)


def output_csv(finite_N, large_N, su3, output_file):
    writer = csv.DictWriter(
        output_file,
        fieldnames=[
            "group_family",
            "Nc",
            "representation",
            "Nf",
            "s0_value",
            "s0_uncertainty",
            "s1_value",
            "s1_uncertainty",
            "s2_value",
            "s2_uncertainty",
        ],
    )
    writer.writeheader()
    for Nc, Nc_data in {**finite_N, **large_N}.items():
        for representation, rep_data in Nc_data.items():
            write_csv_row("Sp", Nc, representation, 0, rep_data, writer)

    write_csv_row("SU", 3, "F", 2, su3, writer)


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
        Nc: compute_s(data) for Nc, data in finite_N_squared_data.items()
    }
    large_N_results = {r"\infty": compute_s(large_N_squared_data)}
    su3_results = {
        r"$SU(3)$, $N_{(\mathrm{f})}=2$ ($m_\pi=139.6\textnormal{ MeV}$)": (
            compute_single_s(su3_squared_data)
        )
    }
    output_table(
        finite_N_results, large_N_results, su3_results, output_file=args.output_table
    )
    if args.output_csv:
        output_csv(
            finite_N_results,
            large_N_results,
            compute_single_s(su3_squared_data),
            output_file=args.output_csv,
        )


if __name__ == "__main__":
    main()
