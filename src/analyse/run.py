import calliope

calliope.set_log_verbosity()


def build_model(model_yaml_path, output_model_path, scenario, run):
    if run is False:
        m = create_inputs(model_yaml_path, scenario)

    m.to_netcdf(output_model_path)
    return None


def create_inputs(model_yaml_path, scenario):
    return calliope.Model(model_yaml_path, scenario=scenario)


if __name__ == '__main__':
    build_model(
        model_yaml_path=snakemake.input.model_yaml_path,
        output_model_path=snakemake.output.output_model_path,
        scenario=snakemake.params.scenario,
        run=snakemake.params.run
    )