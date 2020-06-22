localrules: maps

rule maps:
    message: "Creating Calliope {wildcards.resolution} map"
    input:
        src = "src/analyse/maps.py",
        model = "inputs/{resolution}/model.nc",
        units = eurocalliope("build/data/eurospores/units.geojson")
    output:
        "inputs/{resolution}/map.pdf"
    conda: "../envs/plots.yaml"
    script: "../src/analyse/maps.py"


rule input_netcdf:
    message: "Creating Calliope {wildcards.resolution} map"
    input:
        src = "src/analyse/run.py",
        model_yaml_path = "build/model/{resolution}/model.yaml"
    params:
        scenario = 'directional-rooftop-pv,outer-countries',
        run = False
    output:
        output_model_path = "inputs/{resolution}/model.nc"
    conda: "../envs/calliope.yaml"
    script: "../src/analyse/run.py"
