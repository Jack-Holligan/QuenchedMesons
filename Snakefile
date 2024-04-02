from metadata import metadata

betas = {
    4: [7.62, 7.7, 7.85, 8.0, 8.2],
    6: [15.6, 16.1, 16.5, 16.7, 17.1],
    8: [26.5, 26.7, 26.8, 27.0, 27.3],
}
Ncs = [4, 6, 8]
reps = ["F", "AS", "S"]
channel_names = {
    "PS": "pseudoscalar",
    "S": "scalar",
    "V": "vector",
    "AV": "axialvector",
    "T": "tensor",
    "AT": "axialtensor",
}
channels = list(channel_names.values())
mass_channels = ["vector", "axialvector", "scalar", "tensor", "axialtensor"]
decayconst_channels = ["pseudoscalar", "vector", "axialvector"]
channel_observables = [f"{channel}_masses" for channel in mass_channels] + [
    f"{channel}_decayconsts" for channel in decayconst_channels
]


csvs = expand(
    "csvs/{basename}.csv",
    basename=["w0", "ensemble_masses", "continuum", "large_N", "sumrules"],
)

tables = [
    "assets/tables/w0.tex",
    *expand("assets/tables/chiral_Sp{Nc}.tex", Nc=Ncs),
    "assets/tables/chiral_largeN.tex",
    "assets/tables/sumrules.tex",
]

largeN_plots = expand(
    "assets/plots/largeN_{rep}_{channel_observable}.pdf",
    rep=reps,
    channel_observable=channel_observables,
)
box_plots = expand(
    "assets/plots/{representation}{observable}.pdf",
    representation=["Fundamental", "Antisymmetric", "Symmetric"],
    observable=["Decay", "Mass"],
)
continuum_plots = [
    f"assets/plots/continuum_Sp{Nc}_{channel_observable}_{rep}.pdf"
    for rep in reps
    for channel_observable in channel_observables
    for Nc in Ncs
    if channel_observable.split("_")[0] in metadata.ensembles[Nc][rep]
]
finitesize_plots = expand(
    "assets/plots/Sp6_beta15.6_mF{mass}_finitesize.pdf",
    mass=[-0.8, -0.81],
)

plots = [largeN_plots, box_plots, continuum_plots, finitesize_plots]

definitions = "assets/definitions.tex"


rule all:
    input:
        csvs,
        tables,
        plots,
        definitions,


rule strip_mesons:
    input:
        data="raw_data/Sp{Nc}/beta{beta}/{slug}.out",
        script="src/mesfilter_new.sh",
    output:
        "processed_data/Sp{Nc}/beta{beta}/{slug}/correlators_{slug}.dat",
    shell:
        "cat {input.data} | bash {input.script} > {output}"


rule strip_plaquettes:
    input:
        "raw_data/Sp{Nc}/beta{beta}/{slug}.out",
    output:
        "processed_data/Sp{Nc}/beta{beta}/{slug}/plaquette_{slug}.dat",
    shell:
        """awk '/Plaquette/ {{gsub("Plaquette=", ""); print($8);}}' {input} > {output}"""


rule generate_fit_correlation_function_script:
    input:
        data="src/fit_correlation_function.wls",
        metadata="metadata/ensembles.yaml",
        script="src/wrap_fit_correlation_function.py",
    output:
        "processed_data/Sp{Nc}/beta{beta}/{slug}/fit_correlation_function.wls",
    conda:
        "environment.yml"
    shell:
        "python {input.script} {wildcards.Nc} {wildcards.beta} {wildcards.slug} > {output}"


rule fit_correlation_function:
    input:
        script="processed_data/Sp{Nc}/beta{beta}/{slug}/fit_correlation_function.wls",
        correlators="processed_data/Sp{Nc}/beta{beta}/{slug}/correlators_{slug}.dat",
        plaquettes="processed_data/Sp{Nc}/beta{beta}/{slug}/plaquette_{slug}.dat",
    output:
        expand(
            "processed_data/Sp{{Nc}}/beta{{beta}}/{{slug}}/{{slug}}_{channel}_masses_boots.csv",
            channel=channels,
        ),
        expand(
            "processed_data/Sp{{Nc}}/beta{{beta}}/{{slug}}/{{slug}}_{channel}_decayconsts_boots.csv",
            channel=decayconst_channels,
        ),
        expand(
            "processed_data/Sp{{Nc}}/beta{{beta}}/{{slug}}/{channel}{suffix}.pdf",
            channel=channels,
            suffix=["", "CSD"],
        ),
        expand(
            "processed_data/Sp{{Nc}}/beta{{beta}}/{{slug}}/output_{channel}.txt",
            channel=channels,
        ),
        "processed_data/Sp{Nc}/beta{beta}/{slug}/{slug}.txt",
    log:
        "processed_data/Sp{Nc}/beta{beta}/{slug}/fit_correlation_function.log",
    resources:
        mathematica_licenses=1,
    shell:
        "wolframscript -file {input.script} > {log}"


def ensemble_spectrum_datafiles(wildcards):
    return expand(
        "processed_data/Sp{{Nc}}/beta{{beta}}/{{volume}}B{{beta}}_m{{rep}}{mass}/output_{{channel}}.txt",
        mass=sorted(
            metadata.bare_masses[int(wildcards.Nc)][
                f"{wildcards.volume}B{wildcards.beta}"
            ][wildcards.rep][wildcards.channel],
            reverse=True,
        ),
    )


rule collate_masses:
    input:
        spectrum=ensemble_spectrum_datafiles,
        w0="processed_data/Sp{Nc}/beta{beta}/wflow.dat",
        script="src/collate_masses.py",
    output:
        "processed_data/Sp{Nc}/continuum/{rep}_data/{volume}B{beta}_masses_{channel}_{rep}.txt",
    conda:
        "environment.yml"
    shell:
        "python {input.script} {input.spectrum} --w0_file {input.w0} --output_file {output}"


rule collate_decay_consts:
    input:
        spectrum=ensemble_spectrum_datafiles,
        w0="processed_data/Sp{Nc}/beta{beta}/wflow.dat",
        script="src/collate_masses.py",
    output:
        "processed_data/Sp{Nc}/continuum/{rep}_data/{volume}B{beta}_decayconsts_{channel}_{rep}.txt",
    conda:
        "environment.yml"
    shell:
        "python {input.script} {input.spectrum} --w0_file {input.w0} --output_file {output} --observable decayconst"


rule collate_boots:
    input:
        "processed_data/Sp{Nc}/beta{beta}/{slug}B{beta}_m{rep}{mass}/{slug}B{beta}_m{rep}{mass}_{channel}_{observable}_boots.csv",
    output:
        "processed_data/Sp{Nc}/continuum/{rep}_data/{slug}B{beta}_m{rep}{mass}_{channel}_{observable}_boots.csv",
    shell:
        "cp {input} {output}"


def wflow_log(wildcards):
    ensemble = metadata.flow_ensembles[int(wildcards.Nc)][float(wildcards.beta)]
    return f"raw_data/Sp{{Nc}}/beta{{beta}}/wflow_nc{{Nc}}S{ensemble['Ns']}T{ensemble['Nt']}B{{beta}}.out"


rule wilson_flow:
    input:
        log=wflow_log,
        ensembles="metadata/puregauge.yaml",
        script="src/WilsonFlow.py",
    output:
        text="processed_data/Sp{Nc}/beta{beta}/wflow.dat",
        plot="processed_data/Sp{Nc}/beta{beta}/wflow.pdf",
    log:
        "processed_data/Sp{Nc}/beta{beta}/wflow.log",
    conda:
        "environment.yml"
    shell:
        "python {input.script} --beta {wildcards.beta} --flow_file {input.log} --metadata {input.ensembles} --output_file_main {output.text} --output_file_plot {output.plot} > {log}"


def flow_datafiles(wildcards):
    return expand(
        "processed_data/Sp{{Nc}}/beta{beta}/wflow.dat", beta=betas[int(wildcards.Nc)]
    )


rule combine_wilson_flow:
    input:
        flow_datafiles,
    output:
        "processed_data/Sp{Nc}/wflow.dat",
    shell:
        "cat {input} > {output}"


rule generate_continuum_script:
    input:
        "metadata/continuum.yaml",
        "metadata/puregauge.yaml",
        source="src/continuum.wls",
        script="src/wrap_continuum.py",
    output:
        "processed_data/Sp{Nc}/continuum/{rep}/continuum_{observable}_{channel}.wls",
    conda:
        "environment.yml"
    shell:
        "python {input.script} {wildcards.Nc} {wildcards.rep} {wildcards.channel} {wildcards.observable} > {output}"


def continuum_data(wildcards):
    filelist = []
    for slug in metadata.ensembles[int(wildcards.Nc)][wildcards.rep][wildcards.channel]:
        filelist.append(
            f"processed_data/Sp{{Nc}}/continuum/{{rep}}_data/{slug}_{{observable}}_{{channel}}_{{rep}}.txt"
        )
        filelist.append(
            f"processed_data/Sp{{Nc}}/continuum/{{rep}}_data/{slug}_masses_pseudoscalar_{{rep}}.txt"
        )
        for mass in metadata.bare_masses[int(wildcards.Nc)][slug][wildcards.rep][
            wildcards.channel
        ]:
            filelist.append(
                f"processed_data/Sp{{Nc}}/continuum/{{rep}}_data/{slug}_m{{rep}}{mass}_pseudoscalar_masses_boots.csv"
            )
            filelist.append(
                f"processed_data/Sp{{Nc}}/continuum/{{rep}}_data/{slug}_m{{rep}}{mass}_{{channel}}_{{observable}}_boots.csv"
            )

    return filelist


rule continuum:
    input:
        data=continuum_data,
        script="processed_data/Sp{Nc}/continuum/{rep}/continuum_{observable}_{channel}.wls",
    output:
        "processed_data/Sp{Nc}/continuum/{rep}/{channel}_{observable}_{rep}_Sp{Nc}.pdf",
        "processed_data/Sp{Nc}/continuum/{rep}/{channel}_{observable}_{rep}_Sp{Nc}.dat",
        "processed_data/Sp{Nc}/continuum/{rep}/{channel}_{observable}_{rep}_Sp{Nc}_highbeta.dat",
    log:
        "processed_data/Sp{Nc}/continuum/{rep}/continuum_{observable}_{channel}.log",
    resources:
        mathematica_licenses=1,
    shell:
        "wolframscript -file {input.script} > {log}"


def box_plot_sources(wildcards):
    filelist = []
    for channel in mass_channels:
        for rep in reps:
            if channel in metadata.ensembles[int(wildcards.Nc)][rep]:
                filelist.append(
                    f"processed_data/Sp{{Nc}}/continuum/{rep}/{channel}_masses_{rep}_Sp{{Nc}}.dat"
                )
                if channel in decayconst_channels:
                    filelist.append(
                        f"processed_data/Sp{{Nc}}/continuum/{rep}/{channel}_decayconsts_{rep}_Sp{{Nc}}.dat"
                    )
    return filelist


rule collate_box_plot_inputs:
    input:
        data=box_plot_sources,
        script="src/collate_boxplot_inputs.sh",
    output:
        expand(
            "processed_data/Sp{{Nc}}/continuum/{rep}/{rep}_{observable}.txt",
            observable=["masses", "decayconsts"],
            rep=reps,
        ),
    shell:
        "bash {input.script} processed_data/Sp{wildcards.Nc}/continuum {wildcards.Nc}"


rule box_plot:
    input:
        finite_N_inputs=expand(
            "processed_data/Sp{Nc}/continuum/{rep}/{rep}_{observable}.txt",
            rep=reps,
            observable=["masses", "decayconsts"],
            Nc=Ncs,
        ),
        large_N_inputs=expand(
            "processed_data/largeN/{rep}_{observable}.txt",
            rep=reps,
            observable=["masses", "decayconsts"],
        ),
        script="src/boxplot.wls",
    output:
        boxplots=expand(
            "processed_data/boxplots/{representation}{observable}.pdf",
            representation=["Fundamental", "Antisymmetric", "Symmetric"],
            observable=["Decay", "Mass"],
        ),
        ratioplots=expand(
            "processed_data/boxplots/ftomRatio{rep}.pdf",
            rep=reps,
        ),
    log:
        "processed_data/boxplots/boxplot.log",
    shell:
        "wolframscript -file {input.script} > {log}"


rule collate_box_plots:
    input:
        "processed_data/boxplots/{slug}.pdf",
    output:
        "assets/plots/{slug}.pdf",
    shell:
        "cp {input} {output}"


rule generate_large_N_script:
    input:
        "src/largeN.wls",
    output:
        "processed_data/largeN/largeN_{observable}_{channel}_{rep}.wls",
    shell:
        "sed 's/_SED_REP_/{wildcards.rep}/;s/_SED_CHANNEL_/{wildcards.channel}/;s/_SED_OBSERVABLE_/{wildcards.observable}/' {input} > {output}"


rule large_N:
    input:
        script="processed_data/largeN/largeN_{observable}_{channel}_{rep}.wls",
        data=expand(
            "processed_data/Sp{Nc}/continuum/{{rep}}/{{rep}}_{{observable}}.txt",
            Nc=Ncs,
        ),
    output:
        "processed_data/largeN/{rep}_{channel}_{observable}.pdf",
        "processed_data/largeN/{rep}_{channel}_{observable}.txt",
    log:
        "processed_data/largeN/largeN_{observable}_{channel}_{rep}.log",
    shell:
        "wolframscript -file {input.script} > {log}"


rule collate_masses_large_N:
    input:
        expand(
            "processed_data/largeN/{{rep}}_{channel}_masses.txt",
            channel=mass_channels,
        ),
    output:
        "processed_data/largeN/{rep}_masses.txt",
    shell:
        "awk 1 {input} > {output}"


rule collate_decayconsts_large_N:
    input:
        expand(
            "processed_data/largeN/{{rep}}_{channel}_decayconsts.txt",
            channel=decayconst_channels,
        ),
    output:
        "processed_data/largeN/{rep}_decayconsts.txt",
    shell:
        "awk 1 {input} > {output}"


rule generate_finite_size_script:
    input:
        "src/FiniteSize.wls",
    output:
        "processed_data/Sp{Nc}/beta{beta}/FiniteSize.wls",
    shell:
        "sed 's:_SED_BASEPATH_:processed_data/Sp{wildcards.Nc}/beta{wildcards.beta}:;s:_SED_OUTPUTPATH_:processed_data/Sp{wildcards.Nc}/beta{wildcards.beta}:' {input} > {output}"


rule finite_size:
    input:
        script="processed_data/Sp6/beta15.6/FiniteSize.wls",
        data=expand(
            "processed_data/Sp6/beta15.6/{volume}B15.6_mF{mass}/output_pseudoscalar.txt",
            volume=["S12T24", "S16T32", "S20T40", "S24T48"],
            mass=["-0.8", "-0.81", "-0.82"],
        ),
    output:
        expand(
            "processed_data/Sp6/beta15.6/Sp6_beta15.6_mF{mass}.pdf",
            mass=["-0.8", "-0.81", "-0.82"],
        ),
    log:
        "processed_data/Sp6/beta15.6/finitesize.log",
    shell:
        "wolframscript -file {input.script} > {log}"


rule relocate_finite_size:
    input:
        "processed_data/Sp{Nc}/beta{beta}/Sp{Nc}_beta{beta}_m{Rep}{mass}.pdf",
    output:
        "assets/plots/Sp{Nc}_beta{beta}_m{Rep}{mass}_finitesize.pdf",
    shell:
        "cp {input} {output}"


rule collate_large_N_plots:
    input:
        "processed_data/largeN/{rep}_{channel}_{observable}.pdf",
    output:
        "assets/plots/largeN_{rep}_{channel}_{observable}.pdf",
    shell:
        "cp {input} {output}"


rule collate_contlim_plots:
    input:
        "processed_data/Sp{Nc}/continuum/{rep}/{channel}_{observable}_{rep}_Sp{Nc}.pdf",
    output:
        "assets/plots/continuum_Sp{Nc}_{channel}_{observable}_{rep}.pdf",
    shell:
        "cp {input} {output}"


rule collate_finitesize_plots:
    input:
        "processed_data/Sp{Nc}/beta{beta}/Sp{Nc}_beta{beta}_m{rep}{mass}.pdf",
    output:
        "assets/plots/finitesize_Sp{Nc}_beta{beta}_m{rep}{mass}.pdf",
    shell:
        "cp {input} {output}"


def w0_data(wildcards):
    return [
        f"processed_data/Sp{Nc}/beta{beta}/wflow.dat"
        for Nc, betas in betas.items()
        for beta in betas
    ]


rule w0_table_csv:
    input:
        metadata="metadata/puregauge.yaml",
        data=w0_data,
        script="src/tabulate_w0.py",
    output:
        table="assets/tables/w0.tex",
        csv="csvs/w0.csv",
    conda:
        "environment.yml"
    shell:
        "python {input.script} {input.metadata} --output_table {output.table} --output_csv {output.csv}"


def all_continua(Nc, slug=""):
    return [
        f"processed_data/Sp{Nc}/continuum/{rep}/{channel_observable}_{rep}_Sp{Nc}{slug}.dat"  # fmt: skip (work around Snakefmt bug)
        for rep in reps
        for channel_observable in channel_observables
        if channel_observable.split("_")[0] in metadata.ensembles[Nc][rep]
    ]


def contlim_table_inputs(wildcards):
    return all_continua(int(wildcards.Nc))


rule contlim_tables:
    input:
        data=contlim_table_inputs,
        script="src/LatexChiral.sh",
    output:
        table="assets/tables/chiral_Sp{Nc}.tex",
        caption="processed_data/Sp{Nc}/continuum/caption_vars.tex",
    shell:
        "bash src/LatexChiral.sh {wildcards.Nc} processed_data/Sp{wildcards.Nc}/continuum {output.table} {output.caption}"


rule contlim_csv:
    input:
        data=[all_continua(Nc) for Nc in Ncs]
        + [all_continua(Nc, "_highbeta") for Nc in Ncs],
        script="src/continuum_csv.py",
    output:
        "csvs/continuum.csv",
    conda:
        "environment.yml"
    shell:
        "python {input.script} --output_csv {output}"


def large_N_table_inputs(wildcards):
    return expand(
        "processed_data/largeN/{rep}_{observable}.txt",
        observable=["masses", "decayconsts"],
        rep=reps,
    )


rule large_N_table:
    input:
        data=large_N_table_inputs,
        script="src/LatexChiral_LargeN.sh",
    output:
        table="assets/tables/chiral_largeN.tex",
        caption="processed_data/largeN/caption_vars.tex",
    shell:
        "bash {input.script} processed_data/largeN {output.table} {output.caption}"


rule large_N_csv:
    input:
        data=large_N_table_inputs,
        script="src/largeN_csv.py",
    output:
        "csvs/large_N.csv",
    conda:
        "environment.yml"
    shell:
        "python {input.script} --output_csv {output}"


rule sum_rules:
    input:
        script="src/sumrules.py",
        large_N_data=large_N_table_inputs,
        finite_N_data=[all_continua(Nc) for Nc in Ncs],
    output:
        table="assets/tables/sumrules.tex",
        csv="csvs/sumrules.csv",
    conda:
        "environment.yml"
    shell:
        "python {input.script} --output_table {output.table} --output_csv {output.csv}"


def ensemble_masses(wildcards):
    return [
        "processed_data/Sp{Nc}/beta{beta}/S{NS}T{NT}B{beta}_m{Rep}{m}/output_{channel}.txt".format(
            channel=channel_name, **ensemble
        )
        for ensemble in metadata.metadata.values()
        for channel, channel_name in channel_names.items()
        if f"Use{channel}" in ensemble
    ]


rule ensemble_masses_csv:
    input:
        script="src/ensemble_masses_csv.py",
        metadata="metadata/ensembles.yaml",
        data=ensemble_masses,
    output:
        "csvs/ensemble_masses.csv",
    conda:
        "environment.yml"
    shell:
        "python {input.script} {input.metadata} --output_file {output}"


rule collate_latex_definitions:
    input:
        all_definitions=[
            expand("processed_data/Sp{Nc}/continuum/caption_vars.tex", Nc=Ncs),
            "processed_data/largeN/caption_vars.tex",
        ],
        script="src/collate_latex_definitions.py",
    output:
        definitions,
    conda:
        "environment.yml"
    shell:
        "python {input.script} {input.all_definitions} --output_tex {output}"
