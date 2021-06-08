# coding: utf-8

import argparse
import glob
import os
import re
from pathlib import Path

import calliope
from calliope.core.util.dataset import split_loc_techs
import pandas as pd
import numpy as np
import xarray as xr

ZERO_THRESHOLD = 1e-6  # capacities lower than this will be considered as effectively zero


def save_ds(ds, outfile):
    Path(os.path.dirname(outfile)).mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(outfile)


def combine_spores_to_one_file(indir, new_dimension_name="scenario"):
    input_files = glob.glob(os.path.join(indir, "*.nc"))
    models = {
        spore_search_string(os.path.basename(file)): calliope.read_netcdf(file)
        for file in input_files
    }
    cost_optimal_model = models.pop(0)
    names = cost_optimal_model.inputs.names.to_series()
    input_costs = get_input_costs(cost_optimal_model.inputs)
    output_costs = pd.concat(
        [get_output_cost(model) for model in models.values()],
        keys=models.keys(), names=[new_dimension_name],
    )
    energy_caps = pd.concat(
        [get_caps(model) for model in models.values()],
        keys=models.keys(), names=[new_dimension_name],
    )
    energy_flows = pd.concat(
        [get_flows(model) for model in models.values()],
        keys=models.keys(), names=[new_dimension_name],
    )
    energy_flows = strip_small_techs(energy_caps, energy_flows)
    energy_caps = add_units_to_caps(energy_caps, energy_flows)

    names = names.rename_axis(index="tech").reindex(energy_caps.index.get_level_values("tech").unique())
    ds = (
        xr.Dataset.from_dataframe(energy_caps)
        .merge(xr.Dataset.from_dataframe(energy_flows))
        .merge(xr.Dataset.from_dataframe(input_costs))
    )
    ds["names"] = xr.DataArray.from_series(names)
    ds["cost_result"] = xr.DataArray.from_series(output_costs)

    return ds

def spore_search_string(filename):
    return int(re.search("^spore\\_([\\d]+)\\.nc$", filename).groups()[0])


def get_caps(model):
    return pd.concat([
        clean_series(model._model_data[f"{cap_type}_cap"], keep_demand=False, zero_threshold=ZERO_THRESHOLD)
        .sum(level=["techs", "locs"])
        .rename_axis(index=["tech", "region"])
        .to_frame("{}_capacity".format("nameplate" if cap_type == "energy" else cap_type))
        .div(10)
        for cap_type in ["energy", "storage"]
    ], axis=1)


def get_a_flow(model, flow_direction):
    direction = "out" if flow_direction == "prod" else "in"
    return (
        clean_series(
            model._model_data[f"carrier_{flow_direction}"].sum("timesteps", min_count=1)
        )
        .sum(level=["techs", "carriers", "locs"])
        .rename_axis(index=["tech", "carrier", "region"])
        .div(10)
        .abs()
        .to_frame(f"flow_{direction}")
        .assign(unit='twh_per_year')
        .reset_index()
    )


def get_flows(model):
    prod = get_a_flow(model, "prod")
    con = get_a_flow(model, "con")

    prod.loc[prod.tech.str.find('transport_') > -1, 'unit'] = 'billion_km_per_year'
    con.loc[con.tech.isin(['demand_heavy_transport', 'demand_light_transport']), 'unit'] = 'billion_km_per_year'
    for df in [prod, con]:
        df.loc[df.carrier.isna(), "unit"] = np.nan
        df.loc[df.carrier == 'co2', "unit"] = '100kt_per_year'
    flow = pd.concat(
        [i.set_index(["tech", "region", "carrier", "unit"]).iloc[:, 0] for i in [prod, con]],
        axis=1
    )
    return flow


def clean_series(da, keep_demand=True, zero_threshold=0):
    """
    Get series from dataarray in which we clean out useless info and also sum up all transmission technologies
    """
    initial_series = (
        split_loc_techs(da, return_as="Series")
        .where(lambda x: (~np.isinf(x)) & (abs(x) > zero_threshold))
        .dropna()
    )
    initial_levels = initial_series.index.names
    if "techs" in initial_levels:
        initial_series = split_names_to_clean(initial_series, "techs", ':', initial_levels)
        initial_series = remove_dummy_techs(initial_series, keep_demand)
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


def remove_dummy_techs(series, keep_demand):
    if keep_demand is True:
        series = series.filter(regex="^((?!tech_heat).)*$", axis=0)

    else:
        series = series.filter(regex="^((?!demand).)*$", axis=0)
    heat_storage = series.filter(regex="_heat_storage", axis=0)
    new_tech_names = {
        i: "heat_storage" + i.partition("heat_storage")[-1]
        for i in heat_storage.index.get_level_values('techs')
    }
    _levels = series.index.names
    series = series.append(heat_storage.rename(new_tech_names, level="techs").sum(level=_levels))
    series = series.drop(heat_storage.index.get_level_values("techs").unique(), level="techs")

    if "carriers" in _levels:
        for carrier in ["heat", "transport"]:
            _techs = series.filter(regex="_" + carrier, axis=0)
            series = series.drop(_techs.index)
            series = series.append(_techs.rename(lambda x: x.split("_")[-1], level="carriers").sum(level=_levels))

    return series


def add_units_to_caps(energy_caps, energy_flows):
    multicarrier_info = split_loc_techs(
        cost_optimal_model.inputs.lookup_primary_loc_tech_carriers_out,
        return_as="Series"
    )
    multicarrier_info = (
        multicarrier_info
        .str
        .split("::", expand=True)[2]
        .to_frame("carriers")
    )
    multicarrier_info = (
        multicarrier_info
        .replace({
            x: x.split("_")[1] if "_" in x else x
            for x in np.unique(multicarrier_info.values)
        })
        .groupby(level="techs").first()
        .set_index("carriers", append=True)
        .rename_axis(index=["tech", "carrier"])
    )
    non_multi_carrier = (
        energy_flows.drop(
            multicarrier_info
            .index
            .levels[0],
            level='tech'
        )
        .dropna(subset=["flow_out"])
        .reset_index(["carrier", "unit"])
    )
    levels_to_move = [i for i in energy_flows.index.names if i not in multicarrier_info.index.names]
    multi_carrier = (
        energy_flows
        .unstack(levels_to_move)
        .reindex(multicarrier_info.index)
        .stack(levels_to_move)
        .dropna(subset=["flow_out"])
        .reorder_levels(energy_flows.index.names)
        .reset_index(["carrier", "unit"])
    )
    all_carrier_info = (
        pd.concat([non_multi_carrier, multi_carrier])
        .loc[:, ["carrier", "unit"]]
    )
    all_carrier_info["unit"] = all_carrier_info["unit"].str.replace("per_year", "per_hour")

    energy_caps_with_units = (
        pd.concat([energy_caps, all_carrier_info.reindex(energy_caps.index)], axis=1, sort=True)
        .set_index(["carrier", "unit"], append=True)
    )
    assert len(energy_caps_with_units) == len(energy_caps)

    return energy_caps_with_units


def get_output_cost(model):
    return clean_series(model._model_data.cost.loc[{"costs": "monetary"}])


def get_input_costs(inputs):
    costs = {}
    for var_name, var_data in inputs.data_vars.items():
        if "cost_" not in var_name:
            continue
        if "timesteps" in var_data.dims:
            var_data = var_data.sum("timesteps", min_count=1)
        costs[var_name] = (clean_series(var_data.loc[{"costs": "monetary"}]))
    return pd.concat(costs.values(), keys=costs.keys(), axis=1)


def strip_small_techs(energy_caps, energy_flows):
    levels_to_move = [i for i in energy_flows.index.names if i not in energy_caps.index.names]
    return (
        energy_flows
        .unstack(levels_to_move)
        .loc[energy_caps.index]
        .stack(levels_to_move)
        .reorder_levels(energy_flows.index.names)
        .append(energy_flows.filter(regex="demand", axis=0))
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("indir", help="Input directory in which all datasets to concatenate are found.")
    parser.add_argument("outfile", help="NetCDF file path to save output.")
    parser.add_argument(
        "--scenario_dimension_name", default="scenario",
        help=(
            "Name for the dimension over which all input datasets are concatenated. "
            "If not given, `scenario` will be used."
        )
    )
    args = parser.parse_args()

    ds = combine_spores_to_one_file(args.indir, args.scenario_dimension_name)
    save_ds(ds, args.outfile)
