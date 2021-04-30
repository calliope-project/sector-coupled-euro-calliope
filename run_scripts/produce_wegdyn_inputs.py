# coding: utf-8
import sys
import glob
import os
import re
from pathlib import Path

import calliope
from calliope.core.util.dataset import split_loc_techs
import pandas as pd
import numpy as np
import xarray as xr

ZERO_THRESHOLD = 1e-5  # numbers lower than this will be considered as effectively zero


def combine_spores_to_one_file(indir, outfile, new_dimension_name="spore"):
    input_files = glob.glob(os.path.join(indir, "*.nc"))
    models = {
        spore_search_string(os.path.basename(file)): calliope.read_netcdf(file)
        for file in input_files
    }
    cost_optimal_model = models.pop(0)
    names = cost_optimal_model.inputs.names.to_series()

    energy_caps = pd.concat(
        [get_caps(model) for model in models.values()],
        keys=models.keys(), names=[new_dimension_name],
    )
    energy_flows = pd.concat(
        [get_flows(model) for model in models.values()],
        keys=models.keys(), names=[new_dimension_name],
    )
    names = names.rename_axis(index="tech").reindex(energy_caps.index.get_level_values("tech").unique())
    energy_ds = (
        xr.Dataset.from_dataframe(energy_caps)
        .merge(xr.Dataset.from_dataframe(energy_flows))
    )
    energy_ds["names"] = xr.DataArray.from_series(names)

    Path(os.path.dirname(outfile)).mkdir(parents=True, exist_ok=True)
    energy_ds.to_netcdf(outfile)


def spore_search_string(filename):
    return int(re.search("^spore\\_([\\d]+)\\.nc$", filename).groups()[0])


def get_and_clean_cap(model, cap_type):
    return (
        clean_series(model._model_data[f"{cap_type}_cap"])
        .sum(level=["techs", "locs"])
        .rename_axis(index=["tech", "node"])
        .to_frame(f"{cap_type}_capacity")
        .div(10)
    )


def get_caps(model):
    return pd.concat([
        clean_series(model._model_data[f"{cap_type}_cap"])
        .sum(level=["techs", "locs"])
        .rename_axis(index=["tech", "node"])
        .to_frame(f"{cap_type}_capacity")
        .div(10)
        for cap_type in ["energy", "storage"]
    ], axis=1)


def get_a_flow(model, flow_direction):
    direction = "out" if flow_direction == "prod" else "in"
    return (
        clean_series(model._model_data[f"carrier_{flow_direction}"].sum("timesteps", min_count=1))
        .sum(level=["techs", "carriers", "locs"])
        .rename_axis(index=["tech", "carrier", "node"])
        .div(10)
        .to_frame(f"energy_{direction}")
        .assign(unit='twh_per_year')
        .reset_index()
    )


def get_flows(model):
    prod = get_a_flow(model, "prod")
    con = get_a_flow(model, "con")

    prod.loc[prod.tech.str.find('transport_') > -1, 'unit'] = 'billion_km_per_year'
    con.loc[con.tech.isin(['demand_heavy_transport', 'demand_light_transport']), 'unit'] = 'billion_km_per_year'
    prod = prod.drop(prod[prod.tech.str.find('tech_heat') > -1].index)

    for df in [prod, con]:
        df.loc[df.carrier.isna(), "unit"] = np.nan
        df.loc[df.carrier.str.find('_heat') > -1, f'carrier'] = "heat"
        df.loc[df.carrier == 'co2', "unit"] = '100kt_per_year'
    flow = pd.concat(
        [i.set_index(["tech", "node", "carrier", "unit"]).iloc[:, 0] for i in [prod, con]],
        axis=1
    )
    return flow


def clean_series(da):
    """
    Get series from dataarray in which we clean out useless info and also sum up all transmission technologies
    """
    initial_series = (
        split_loc_techs(da, return_as="Series")
        .where(lambda x: (~np.isinf(x)) & (abs(x) > ZERO_THRESHOLD))
        .dropna()
    )
    initial_levels = initial_series.index.names
    if "techs" in initial_levels:
        initial_series = split_names_to_clean(initial_series, "techs", ':', initial_levels)
    if "locs" in initial_levels:
        initial_series = split_names_to_clean(initial_series, "locs", '_', initial_levels)

    return initial_series


def split_names_to_clean(series, level_of_interest, split_char, initial_levels):
    non_levels = [i for i in series.index.names if i != level_of_interest]
    focussed_series = series.unstack(non_levels)
    focussed_series.index = focussed_series.index.str.split(split_char, expand=True)
    return (
        focussed_series
        .groupby(level=0).sum(min_count=1)
        .rename_axis(index=level_of_interest)
        .unstack()
        .reorder_levels(initial_levels)
        .dropna()
    )


if __name__ == "__main__":

    combine_spores_to_one_file(*sys.argv[1:])
