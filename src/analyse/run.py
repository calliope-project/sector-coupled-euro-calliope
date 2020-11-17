import calliope

calliope.set_log_verbosity()


def build_model(model_yaml_path, output_model_path, scenario, run):
    if run is False:
        m = create_inputs(model_yaml_path, scenario)
        # TODO: fix this in calliope, rather than here
        m._model_data["group_carrier_prod_max"] = m._model_data["group_carrier_prod_max"].astype(float)
        m._model_data["group_carrier_prod_max"].attrs['is_result'] = 0

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