channels = ["pseudoscalar", "vector", "axialvector", "scalar", "tensor", "axialtensor"]


rule all:
    input:
        # "processed_data/Sp4/beta7.62/S24T48B7.62_mAS-1.13/S24T48B7.62_mAS-1.13.txt"
        "processed_data/Sp4/beta7.62/wflow.pdf"

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


rule WilsonFlow:
    input:
        log = "raw_data/Sp{Nc}/beta{beta}/out_wflow",
        ensembles = "metadata/puregauge.yaml"
    output:
        text = "processed_data/Sp{Nc}/beta{beta}/wflow.dat",
        plot = "processed_data/Sp{Nc}/beta{beta}/wflow.pdf"
    shell:
        "python src/WilsonFlow.py --beta {wildcards.beta} --flow_file {input.log} --metadata {input.ensembles} --output_file_main {output.text} --output_file_plot {output.plot}"
