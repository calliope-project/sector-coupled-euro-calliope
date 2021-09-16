import os

import pandas as pd
import calliope

import visualisation_util as util


def consolidate_spores(
    path_to_cost_optimal_model, path_to_spore,
    technical_potential_area, technical_potential_protected_area,
    output_dir, config
):
    potential_area = {
        "technical-potential": pd.read_csv(technical_potential_area, index_col=0, header=0),
        "technical-potential-protected": pd.read_csv(technical_potential_protected_area, index_col=0, header=0)
    }
    cost_opt_model = calliope.read_netcdf(path_to_cost_optimal_model)
    spore_model = calliope.read_netcdf(path_to_spore)
    print("Loaded models")

    spore_util = util.VisUtil(
        "spore", spore_model, config, potential_area, inputs=cost_opt_model.inputs
    )

    _spores_to_csv(spore_util, output_dir)


def _spores_to_csv(util, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    metrics = {
        k: v if isinstance(v, pd.Series) else v.stack()
        for k, v in util.metrics.items()
        if (isinstance(v, pd.DataFrame) or isinstance(v, pd.Series))
    }

    for k, v in metrics.items():
        v.rename(k).to_csv(os.path.join(output_dir, f"{k}.csv"))


if __name__ == '__main__':
    consolidate_spores(
        path_to_cost_optimal_model=snakemake.input.cost_opt_model,
        path_to_spore=snakemake.input.spore,
        technical_potential_area=snakemake.input.technical_potential_area,
        technical_potential_protected_area=snakemake.input.technical_potential_protected_area,
        config=snakemake.params.config,
        output_dir=snakemake.output[0]
    )
