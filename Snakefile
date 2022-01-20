
PYTHON_SCRIPT = "PYTHONPATH=./ python {input} {output}"
PYTHON_SCRIPT_WITH_CONFIG = PYTHON_SCRIPT + " {CONFIG_FILE}"

configfile: "config/default.yaml"
include: "rules/construct.smk"
include: "rules/analyse.smk"
include: "rules/sync.smk"
include: "rules/run.smk"

localrules: all, clean, make_runs
onstart:
    shell("mkdir -p build/logs build/eurospores build/national build/model")

onsuccess:
    if "email" in config.keys():
        shell("echo "" | mail -s 'eurospores succeeded' {config[email]}")
onerror:
    if "email" in config.keys():
        shell("echo "" | mail -s 'eurospores crashed' {config[email]}")
wildcard_constraints:
    resolution = "((national)|(eurospores))", # supported spatial resolutions


rule all:
    message: "Prepare EuroSPORES model runs."
    input:
        "run_scripts/eurospores.sh",
        "build/figures/eurospores/map.pdf"


rule generate_pre_builds:
    message: "Generate pre-built models for all weather years"
    input:
        expand(
            "build/model/{resolution}/model-{year}.yaml",
            resolution=["national", "eurospores"], year=[i for i in range(2010, 2019)]
        )
    output: "build/pre-built-model-{date}.zip"
    shell:
        """
        pushd build
        zip -r $OLDPWD/{output} model/ */annual-demand* annual_industry_energy_demand*  -x "*storyline.yaml" "*bau*" "*water-heat*" "*gshp*" "*ashp*" "*space-heat*"
        popd
        """


rule make_runs:
    message: "Creating Calliope {wildcards.resolution} run scripts"
    input:
        model = "build/model/{resolution}/model.yaml"
    output:
        run_script = "run_scripts/{resolution}.sh"
    conda: "envs/calliope.yaml"
    shell: "calliope generate_runs --kind bsub --scenarios directional-rooftop-pv --cluster_threads 4 --cluster_mem 50G --cluster_time 240 ../{input.model} {output.run_script}"


rule clean: # removes all generated results
    shell:
        """
        rm -r ./build/*
        echo "Data downloaded to data/ has not been cleaned."
        """

