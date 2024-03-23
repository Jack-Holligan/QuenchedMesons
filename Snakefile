from metadata import metadata

betas = {
    4: [7.62, 7.7, 7.85, 8.0, 8.2],
    6: [15.6, 16.1, 16.5, 16.7, 17.1],
    8: [26.5, 26.7, 26.8, 27.0, 27.3],
}

channels = ["pseudoscalar", "vector", "axialvector", "scalar", "tensor", "axialtensor"]


rule all:
    input:
        # "processed_data/Sp4/beta7.62/S24T48B7.62_mAS-1.13/S24T48B7.62_mAS-1.13.txt"
        # "processed_data/Sp4/beta7.62/wflow.pdf"
        # "processed_data/Sp4/wflow.dat"
        # "processed_data/Sp4/continuum/AS_data/S24T48B7.62_massPS_AS.txt"
        "processed_data/Sp4/continuum/F/vector_mass_F_Sp4.dat"

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
    shell:
        "python src/wrap_fit_correlation_function.py {wildcards.Nc} {wildcards.beta} {wildcards.slug} > {output}"

rule fit_correlation_function:
    input:
        script = "processed_data/Sp{Nc}/beta{beta}/{slug}/fit_correlation_function.wls",
        correlators = "processed_data/Sp{Nc}/beta{beta}/{slug}/correlators_{slug}.dat",
        plaquettes = "processed_data/Sp{Nc}/beta{beta}/{slug}/plaquette_{slug}.dat"
    output:
        expand(
            "processed_data/Sp{{Nc}}/beta{{beta}}/{{slug}}/{{slug}}_{channel}_mass_boots.csv",
            channel=channels,
        ),
        expand(
            "processed_data/Sp{{Nc}}/beta{{beta}}/{{slug}}/{{slug}}_{channel}_decayconst_boots.csv",
            channel=channels[:3],
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
        "processed_data/Sp{Nc}/beta{beta}/{slug}/fit_correlation_function.wls"
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
        "processed_data/Sp{Nc}/continuum/{rep}_data/{volume}B{beta}_mass_{channel}_{rep}.txt"
    shell:
        "python src/collate_masses.py {input.spectrum} --w0_file {input.w0} --output_file {output}"

rule CollateDecayConsts:
    input:
        spectrum = ensemble_spectrum_datafiles,
        w0 = "processed_data/Sp{Nc}/beta{beta}/wflow.dat"
    output:
        "processed_data/Sp{Nc}/continuum/{rep}_data/{volume}B{beta}_decayconst_{channel}_{rep}.txt"
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
    shell:
        "python src/wrap_continuum.py {wildcards.Nc} {wildcards.rep} {wildcards.channel} {wildcards.observable} > {output}"

def continuum_data(wildcards):
    filelist = []
    for slug in metadata.ensembles[int(wildcards.Nc)][wildcards.rep][wildcards.channel]:
        filelist.append(f"processed_data/Sp{{Nc}}/continuum/{{rep}}_data/{slug}_{{observable}}_{{channel}}_{{rep}}.txt")
        filelist.append(f"processed_data/Sp{{Nc}}/continuum/{{rep}}_data/{slug}_mass_pseudoscalar_{{rep}}.txt")
        for mass in metadata.bare_masses[int(wildcards.Nc)][slug][wildcards.rep][wildcards.channel]:
            filelist.append(f"processed_data/Sp{{Nc}}/continuum/{{rep}}_data/{slug}_m{{rep}}{mass}_pseudoscalar_mass_boots.csv")
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
        "processed_data/Sp{Nc}/continuum/{rep}/continuum_{observable}_{channel}.wls"
    resources:
        mathematica_licenses = 1
    shell:
        "wolframscript -file {input.script} > {log}"

rule CollateBoxplotInputs:
    input:
        expand(
            "processed_data/Sp{{Nc}}/continuum/{{rep}}/mass_{channel}_{{rep}}_Sp{{Nc}}.dat",
            channel=["vector axialvector scalar tensor axialtensor"],
        ),
        expand(
            "processed_data/Sp{{Nc}}/continuum/{{rep}}/decayconst_{channel}_{{rep}}_Sp{{Nc}}.dat",
            channel=["pseudoscalar vector axialvector"],
        )
    output:
        expand(
            "processed_data/Sp{{Nc}}/continuum/{{rep}}/{{rep}}_{observable}.txt",
            observable=["masses", "decayconsts"],
        )
    shell:
        "bash src/collate_boxplot_inputs.sh processed_data/Sp{Nc}/continuum {Nc}"

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
            "processed_data/Sp{{Nc}}/{rep}/{rep}_{observable}.txt",
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
