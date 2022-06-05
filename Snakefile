from datetime import date
import os

configfile: "config/default.yaml"
configfile: f"config/{config['projection_year']}.yaml"

root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
script_dir = f"{root_dir}src/"
template_dir = f"{root_dir}src/template/"

module eurocalliope:
    snakefile: config["euro-calliope-snakefile"]
    config: config["euro-calliope"]

use rule * from eurocalliope as ec_*

include: "rules/construct.smk"
include: "rules/analyse.smk"
include: "rules/sync.smk"
include: "rules/run.smk"

localrules: all, clean
localrules: ec_download_eez, ec_download_stations_database, ec_download_raw_nuts_units, ec_download_raw_load, ec_download_basins_database, ec_download_raw_gadm_administrative_borders, ec_download_potentials, ec_stations_database, ec_raw_gadm_administrative_borders, ec_basins_database, ec_download_capacity_factors_wind_and_solar

onstart:
    shell("mkdir -p build/logs build/data/ehighways build/data/national build/model")

onsuccess:
    if "email" in config.keys():
        shell("echo "" | mail -s 'sector-coupled euro-calliope succeeded' {config[email]}")
onerror:
    if "email" in config.keys():
        shell("echo "" | mail -s 'sector-coupled euro-calliope crashed' {config[email]}")


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
