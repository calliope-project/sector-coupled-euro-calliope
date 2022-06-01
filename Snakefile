from datetime import date

configfile: "config/default.yaml"
euro_calliope_config = f"./config/euro-calliope-{config['projection_year']}.yaml"
configfile: euro_calliope_config

subworkflow eurocalliope:
    workdir: "./euro-calliope"
    snakefile: "./euro-calliope/Snakefile"
    configfile: euro_calliope_config

include: "rules/construct.smk"
include: "rules/analyse.smk"
include: "rules/sync.smk"
include: "rules/run.smk"

localrules: all, clean
onstart:
    shell("mkdir -p build/logs build/ehighways build/national build/model")

onsuccess:
    if "email" in config.keys():
        shell("echo "" | mail -s 'sector-coupled euro-calliope succeeded' {config[email]}")
onerror:
    if "email" in config.keys():
        shell("echo "" | mail -s 'sector-coupled euro-calliope crashed' {config[email]}")
wildcard_constraints:
    resolution = "((national)|(ehighways))", # supported spatial resolutions


rule all:
    message: "Prepare sector-coupled Euro-Calliope model runs at e-highways and national spatial resolution."
    input:
        "build/figures/ehighways/map.pdf",
        f"build/pre-built-model-{date.today()}.zip"


rule generate_pre_builds:
    message: "Generate pre-built models for all weather years"
    input:
        expand(
            "build/model/{resolution}/model-{year}.yaml",
            resolution=["national", "ehighways"], year=[i for i in range(2010, 2019)]
        )
    output: "build/pre-built-model-{date}.zip"
    shell:
        """
        pushd build
        zip -r $OLDPWD/{output} model/ */annual-demand* annual_industry_energy_demand*  -x "*bau*" "*water-heat*" "*gshp*" "*ashp*" "*space-heat*"
        popd
        """


rule clean: # removes all generated results
    shell:
        """
        rm -r ./build/*
        echo "Data downloaded to data/ has not been cleaned."
        """


rule clean_euro_calliope: # removes all generated results in the subworkflow
    shell:
        """
        rm -r ./euro-calliope/build/*
        echo "Data downloaded to euro-calliope/data/ has not been cleaned."
        """
