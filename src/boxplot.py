#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib import ticker

from numpy import nan
import pandas as pd


shortnames = {
    "pseudoscalar": "PS",
    "vector": "V",
    "axialvector": "AV",
    "tensor": "T",
    "axialtensor": "AT",
    "scalar": "S",
    "mass": "m",
    "decayconst": "f",
}

channels = {
    "mass": ["vector", "axialvector", "scalar", "tensor", "axialtensor"],
    "decayconst": ["pseudoscalar", "vector", "axialvector"],
}

colours = {
    "4": "purple",
    "6": "green",
    "8": "blue",
    r"\infty": "darkorange",
}

channel_offsets = {
    "mass_hat_squared": {
        "vector": 0,
        "axialvector": 1,
        "scalar": 2,
        "tensor": 3,
        "axialtensor": 4,
    },
    "decayconst_hat_squared_over_Nc": {
        "pseudoscalar": 0,
        "vector": 1,
        "axialvector": 2,
    },
}


def get_data():
    all_finite_N = pd.read_csv("csvs/continuum.csv")
    finite_N = all_finite_N[all_finite_N.note == "all beta"]

    large_N = pd.read_csv("csvs/large_N.csv")
    large_N["Nc"] = r"\infty"
    large_N["chiral_limit_value"] = large_N.large_N_value
    large_N["chiral_limit_uncertainty"] = large_N.large_N_uncertainty

    return pd.concat([finite_N, large_N], ignore_index=True)


def get_label(channel, representation):
    shortname = shortnames[channel]
    if representation == "fundamental":
        return f"$\\mathrm{{{shortname}}}$"
    elif representation == "antisymmetric":
        return f"$\\mathrm{{{shortname.lower()}}}$"
    elif representation == "symmetric":
        return f"$\\mathcal{{{shortname}}}$"
    else:
        raise ValueError(f"Don't know what {represenation=} means.")


def offset(datum):
    Nc_offsets = {"4": 0.2, "6": 0.4, "8": 0.6, r"\infty": 0.8}
    return channel_offsets[datum.observable][datum.channel] + Nc_offsets[str(datum.Nc)]


def add_ticks(ax, observable, representation):
    for gridline_index in range(1, len(channels[observable])):
        ax.axvline(gridline_index, alpha=0.2, color="black")

    xtick_positions = [
        0.5 + tick_index for tick_index, _ in enumerate(channels[observable])
    ]
    xtick_labels = [
        get_label(channel, representation) for channel in channels[observable]
    ]
    ax.set_xticks(xtick_positions, labels=xtick_labels)

    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.grid(which="both", axis="y", dashes=(1, 4), lw=0.5, alpha=0.5)


def plot_single(data, filename, representation, observable):
    fig, ax = plt.subplots(figsize=(5, 3), layout="constrained")

    symbol = f"$\\hat{{{shortnames[observable]}}}_{{\\chi}}$"
    ax.set_ylabel(symbol, rotation=0, labelpad=10)

    bar_thickness = 40 / len(channels[observable])
    style = {
        "elinewidth": bar_thickness,
        "ls": "none",
        "capsize": 1.1 * bar_thickness / 2,
        "capthick": 1,
    }

    for Nc in "4", "6", "8", r"\infty":
        ax.errorbar(
            [nan],
            [nan],
            yerr=[nan],
            color=colours[Nc],
            label=f"$\\mathrm{{Sp}}({Nc})$",
            **style,
        )

    for _, datum in data.iterrows():
        ax.errorbar(
            [offset(datum)],
            [datum.chiral_limit_value],
            yerr=[datum.chiral_limit_uncertainty],
            color=colours[str(datum.Nc)],
            **style,
        )

    add_ticks(ax, observable, representation)

    ax.axhline(0, color="black", lw=0.6)
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)

    fig.savefig(filename)
    plt.close(fig)


def main():
    plt.style.use("styles/paper.mplstyle")
    data = get_data()
    for observable in "mass", "decayconst":
        for representation in "fundamental", "antisymmetric", "symmetric":
            subset_data = data[
                (data.representation == representation)
                & (data.observable.str.startswith(observable))
            ]
            filename = f"processed_data/boxplots/{representation}_{observable}.pdf"
            plot_single(subset_data, filename, representation, observable)


if __name__ == "__main__":
    main()
