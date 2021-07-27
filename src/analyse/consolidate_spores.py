import glob
import os

import pandas as pd

import calliope
from friendly_calliope.io import write_dpkg

import visualisation_util as util


def consolidate_spores(
    spores_dir, technical_potential_area, technical_potential_protected_area,
    output_dir, scenario, config
):
    scenario_utils = []
    cost_opt_models = glob.glob(os.path.join(spores_dir, "**", "spore_0.nc"))
    cost_opt_model = calliope.read_netcdf(cost_opt_models[0])
    print("loaded cost optimal model")
    potential_area = {
        "technical-potential": pd.read_csv(technical_potential_area, index_col=0, header=0),
        "technical-potential-protected": pd.read_csv(technical_potential_protected_area, index_col=0, header=0)
    }
    cost_opt_utils = util.VisUtil('cost_opt', cost_opt_model, config, potential_area)
    model_files = [
        i for i in
        glob.glob(os.path.join(spores_dir, f"*{scenario}*", "spore_*.nc"))
        if "spore_0.nc" not in i
    ]
    processed_spores = {}
    spore_num = 0
    for file in model_files:
        scenario_utils.append(
            util.VisUtil(
                spore_num, calliope.read_netcdf(file), config, potential_area,
                inputs=cost_opt_model.inputs
            )
        )
        processed_spores[spore_num] = file
        spore_num += 1
        if spore_num % 20 == 0:
            print(f"{spore_num / len(model_files) * 100:.2f}% through consolidating spores")

    grouped_metrics = util.get_grouped_metrics(scenario_utils, "scenario")
    print("loaded all spores metrics")

    meta_spores = {
        "name": "euro-calliope-spores-results",
        "description": "Calliope Euro-SPORES output dataset",
        "keywords": ["calliope", "SPORES", "2h_resolution"],
        "license": "CC-BY-4.0"
    }
    spores_metrics = {
        k: v if isinstance(v, pd.Series) else v.stack()
        for k, v in grouped_metrics.items()
        if (isinstance(v, pd.DataFrame) or isinstance(v, pd.Series))
    }
    write_dpkg(
        spores_metrics, os.path.join(output_dir, "spores"), meta_spores
    )
    print("saved spores metrics")

    meta_cost_opt = {
        "name": "euro-calliope-cost-opt-data",
        "description": "Calliope SPORES cost optimal baseline output dataset",
        "keywords": ["calliope", "SPORES", "2h_resolution", "cost-optimal"],
        "license": "CC-BY-4.0"
    }
    cost_opt_metrics = {
        k: v if isinstance(v, pd.Series) else v.stack()
        for k, v in cost_opt_utils.metrics.items()
        if (isinstance(v, pd.DataFrame) or isinstance(v, pd.Series))
    }
    # Hack to add an additional dimension here so it keeps `friendly_data` happy
    cost_opt_metrics["transmission_flows"] = (
        cost_opt_metrics["transmission_flows"]
        .to_frame("cost_opt")
        .rename_axis(columns="scenario")
        .stack()
    )
    write_dpkg(
        cost_opt_metrics, os.path.join(output_dir, "cost_opt"), meta_cost_opt
    )
    print("saved cost optimal metrics")

    pd.Series(processed_spores).to_csv(os.path.join(output_dir, "processed_spores.csv"))


if __name__ == '__main__':
    consolidate_spores(
        spores_dir=snakemake.input.spores_dir,
        technical_potential_area=snakemake.input.technical_potential_area,
        technical_potential_protected_area=snakemake.input.technical_potential_protected_area,
        output_dir=snakemake.output[0],
        config=snakemake.params.config,
        scenario=snakemake.wildcards.scenario
    )
