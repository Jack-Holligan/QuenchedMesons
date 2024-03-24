from metadata import metadata

betas = {
    4: [7.62, 7.7, 7.85, 8.0, 8.2],
    6: [15.6, 16.1, 16.5, 16.7, 17.1],
    8: [26.5, 26.7, 26.8, 27.0, 27.3],
}
Ncs = [4, 6, 8]
reps = ["F", "AS", "S"]
channels = ["pseudoscalar", "vector", "axialvector", "scalar", "tensor", "axialtensor"]
mass_channels = ["vector", "axialvector", "scalar", "tensor", "axialtensor"]
decayconst_channels = ["pseudoscalar", "vector", "axialvector"]
channel_observables = [
    f"{channel}_masses" for channel in mass_channels
] + [
    f"{channel}_decayconsts" for channel in decayconst_channels
]

rule all:
    input:
        expand(
            "processed_data/largeN/{rep}_{channel_observable}.pdf",
            rep=reps,
            channel_observable=channel_observables,
        ),
        expand(
            "processed_data/Sp{Nc}/continuum/chiral_mass_Sp{Nc}.pdf",
            Nc=Ncs,
        ),
        expand(
            "processed_data/Sp{Nc}/continuum/chiral_decayconst_Sp{Nc}.pdf",
            Nc=Ncs,
        ),
        [
            f"processed_data/Sp{Nc}/continuum/{rep}/{channel_observable}_{rep}_Sp{Nc}.pdf"
            for rep in reps
            for channel_observable in channel_observables
            for Nc in Ncs
            if channel_observable.split("_")[0] in metadata.ensembles[Nc][rep]
        ]

rule strip_mesons:
    input:
        "raw_data/Sp{Nc}/beta{beta}/{slug}.out"
    output:
        "processed_data/Sp{Nc}/beta{beta}/{slug}/correlators_{slug}.dat"
    shell:
        "cat {input} | ./src/mesfilter_new.sh > {output}"

rule strip_plaquettes:
    input:
        "raw_data/Sp{Nc}/beta{beta}/{slug}.out"
    output:
        "processed_data/Sp{Nc}/beta{beta}/{slug}/plaquette_{slug}.dat"
    shell:
        """awk '/Plaquette/ {{gsub("Plaquette=", ""); print($8);}}' {input} > {output}"""

rule generate_fit_correlation_function_script:
    input:
        "src/fit_correlation_function.wls",
        "metadata/ensembles.yaml"
    output:
        "processed_data/Sp{Nc}/beta{beta}/{slug}/fit_correlation_function.wls"
    conda:
        "environment.yml"
    shell:
        "python src/wrap_fit_correlation_function.py {wildcards.Nc} {wildcards.beta} {wildcards.slug} > {output}"

rule fit_correlation_function:
    input:
        script = "processed_data/Sp{Nc}/beta{beta}/{slug}/fit_correlation_function.wls",
        correlators = "processed_data/Sp{Nc}/beta{beta}/{slug}/correlators_{slug}.dat",
        plaquettes = "processed_data/Sp{Nc}/beta{beta}/{slug}/plaquette_{slug}.dat"
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
        "processed_data/Sp{Nc}/beta{beta}/{slug}/{slug}.txt"
    log:
        "processed_data/Sp{Nc}/beta{beta}/{slug}/fit_correlation_function.log"
    resources:
        mathematica_licenses = 1
    shell:
        "wolframscript -file {input.script} > {log}"


def ensemble_spectrum_datafiles(wildcards):
    return expand(
        "processed_data/Sp{{Nc}}/beta{{beta}}/{{volume}}B{{beta}}_m{{rep}}{mass}/output_{{channel}}.txt",
        mass=metadata.bare_masses[int(wildcards.Nc)][f"{wildcards.volume}B{wildcards.beta}"][wildcards.rep][wildcards.channel],
    )

rule CollateMasses:
    input:
        spectrum = ensemble_spectrum_datafiles,
        w0 = "processed_data/Sp{Nc}/beta{beta}/wflow.dat"
    output:
        "processed_data/Sp{Nc}/continuum/{rep}_data/{volume}B{beta}_masses_{channel}_{rep}.txt"
    conda:
        "environment.yml"
    shell:
        "python src/collate_masses.py {input.spectrum} --w0_file {input.w0} --output_file {output}"

rule CollateDecayConsts:
    input:
        spectrum = ensemble_spectrum_datafiles,
        w0 = "processed_data/Sp{Nc}/beta{beta}/wflow.dat"
    output:
        "processed_data/Sp{Nc}/continuum/{rep}_data/{volume}B{beta}_decayconsts_{channel}_{rep}.txt"
    conda:
        "environment.yml"
    shell:
        "python src/collate_masses.py {input.spectrum} --w0_file {input.w0} --output_file {output} --observable decayconst"

rule CollateBoots:
    input:
        "processed_data/Sp{Nc}/beta{beta}/{slug}B{beta}_m{rep}{mass}/{slug}B{beta}_m{rep}{mass}_{channel}_{observable}_boots.csv"
    output:
        "processed_data/Sp{Nc}/continuum/{rep}_data/{slug}B{beta}_m{rep}{mass}_{channel}_{observable}_boots.csv"
    shell:
        "cp {input} {output}"

rule WilsonFlow:
    input:
        log = "raw_data/Sp{Nc}/beta{beta}/out_wflow",
        ensembles = "metadata/puregauge.yaml"
    output:
        text = "processed_data/Sp{Nc}/beta{beta}/wflow.dat",
        plot = "processed_data/Sp{Nc}/beta{beta}/wflow.pdf"
    conda:
        "environment.yml"
    shell:
        "python src/WilsonFlow.py --beta {wildcards.beta} --flow_file {input.log} --metadata {input.ensembles} --output_file_main {output.text} --output_file_plot {output.plot}"


def flow_datafiles(wildcards):
    return expand(
        "processed_data/Sp{{Nc}}/beta{beta}/wflow.dat",
        beta=betas[int(wildcards.Nc)]
    )

rule CombineWilsonFlow:
    input:
        flow_datafiles
    output:
        "processed_data/Sp{Nc}/wflow.dat"
    shell:
        "cat {input} > {output}"

rule GenerateContinuumScript:
    input:
        "metadata/ensembles.yaml",
        "metadata/puregauge.yaml",
        "src/continuum.wls"
    output:
        "processed_data/Sp{Nc}/continuum/{rep}/continuum_{observable}_{channel}.wls"
    conda:
        "environment.yml"
    shell:
        "python src/wrap_continuum.py {wildcards.Nc} {wildcards.rep} {wildcards.channel} {wildcards.observable} > {output}"

def continuum_data(wildcards):
    filelist = []
    for slug in metadata.ensembles[int(wildcards.Nc)][wildcards.rep][wildcards.channel]:
        filelist.append(f"processed_data/Sp{{Nc}}/continuum/{{rep}}_data/{slug}_{{observable}}_{{channel}}_{{rep}}.txt")
        filelist.append(f"processed_data/Sp{{Nc}}/continuum/{{rep}}_data/{slug}_masses_pseudoscalar_{{rep}}.txt")
        for mass in metadata.bare_masses[int(wildcards.Nc)][slug][wildcards.rep][wildcards.channel]:
            filelist.append(f"processed_data/Sp{{Nc}}/continuum/{{rep}}_data/{slug}_m{{rep}}{mass}_pseudoscalar_masses_boots.csv")
            filelist.append(f"processed_data/Sp{{Nc}}/continuum/{{rep}}_data/{slug}_m{{rep}}{mass}_{{channel}}_{{observable}}_boots.csv")

    return filelist

rule Continuum:
    input:
        data = continuum_data,
        script = "processed_data/Sp{Nc}/continuum/{rep}/continuum_{observable}_{channel}.wls"
    output:
        "processed_data/Sp{Nc}/continuum/{rep}/{channel}_{observable}_{rep}_Sp{Nc}.pdf",
        "processed_data/Sp{Nc}/continuum/{rep}/{channel}_{observable}_{rep}_Sp{Nc}.dat"
    log:
        "processed_data/Sp{Nc}/continuum/{rep}/continuum_{observable}_{channel}.log"
    resources:
        mathematica_licenses = 1
    shell:
        "wolframscript -file {input.script} > {log}"

def box_plot_sources(wildcards):
    filelist = []
    for channel in mass_channels:
        for rep in reps:
            if channel in metadata.ensembles[int(wildcards.Nc)][rep]:
               filelist.append(f"processed_data/Sp{{Nc}}/continuum/{rep}/{channel}_masses_{rep}_Sp{{Nc}}.dat")
               if channel in decayconst_channels:
                  filelist.append(f"processed_data/Sp{{Nc}}/continuum/{rep}/{channel}_decayconsts_{rep}_Sp{{Nc}}.dat")
    return filelist

rule CollateBoxplotInputs:
    input:
        box_plot_sources
    output:
        expand(
            "processed_data/Sp{{Nc}}/continuum/{rep}/{rep}_{observable}.txt",
            observable=["masses", "decayconsts"],
            rep=reps,
        )
    shell:
        "bash src/collate_boxplot_inputs.sh processed_data/Sp{wildcards.Nc}/continuum {wildcards.Nc}"

rule GenerateBoxplotScript:
    input:
        "src/boxplot.wls"
    output:
        "processed_data/Sp{Nc}/continuum/boxplot.wls"
    shell:
        "sed 's/_SED_NC_/{wildcards.Nc}/' {input} > {output}"

rule Boxplot:
    input:
        inputs = expand(
            "processed_data/Sp{{Nc}}/continuum/{rep}/{rep}_{observable}.txt",
            rep=["F", "AS", "S"],
            observable=["masses", "decayconsts"],
        ),
        script = "processed_data/Sp{Nc}/continuum/boxplot.wls"
    output:
        "processed_data/Sp{Nc}/continuum/chiral_mass_Sp{Nc}.pdf",
        "processed_data/Sp{Nc}/continuum/chiral_decayconst_Sp{Nc}.pdf"
    log:
        "processed_data/Sp{Nc}/continuum/boxplot.log"
    shell:
        "wolframscript -file {script} > {log}"

rule GenerateLargeNScript:
    input:
        "src/largeN.wls"
    output:
        "processed_data/largeN/largeN_{observable}_{channel}_{rep}.wls"
    shell:
        "sed 's/_SED_REP_/{wildcards.rep}/;s/_SED_CHANNEL_/{wildcards.channel}/;s/_SED_OBSERVABLE_/{wildcards.observable}/' {input} > {output}"

rule LargeN:
    input:
        "processed_data/largeN/largeN_{observable}_{channel}_{rep}.wls",
        expand(
            "processed_data/Sp{Nc}/Continuum/{{rep}}/{{rep}}_{{observable}}.txt"
        )
    output:
        "processed_data/largeN/{rep}_{channel}_{observable}.pdf",
        "processed_data/largeN/{rep}_{channel}_{observable}.txt"
    log:
        "processed_data/largeN/largeN_{observable}_{channel}_{rep}.log"
    rule:
        "wolframscript -file {script} > {log}"

rule CollateMassesLargeN:
    input:
        expand(
            "processed_data/largeN/{{rep}}_{channel}_masses.txt",
            channel=mass_channels,
        )
    output:
        "processed_data/largeN/{rep}_masses.txt"
    rule:
        "cat {input} > {output}"

rule CollateDecayConstsLargeN:
    input:
        expand(
            "processed_data/largeN/{{rep}}_{channel}_decayconsts.txt",
            channel=decayconst_channels,
        )
    output:
        "processed_data/largeN/{rep}_decayconsts.txt"
    shell:
        "cat {input} > {output}"
