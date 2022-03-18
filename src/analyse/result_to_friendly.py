import calliope

from friendly_calliope.io import write_dpkg
from friendly_calliope import consolidate_calliope_output

SPORE = 1

import sys
import os

import pandas as pd

sys.path.append("src/construct")
import util

ENERGY_PRODUCERS = {
    "waste_supply": "Waste",
    "biofuel_supply": "Biofuels",
    "hydro_reservoir": "Renewable electricity",
    "hydro_run_of_river": "Renewable electricity",
    "nuclear": "Nuclear electricity",
    "open_field_pv": "Renewable electricity",
    "roof_mounted_pv": "Renewable electricity",
    "wind_offshore": "Renewable electricity",
    "wind_onshore": "Renewable electricity",
}

HEAT_TECHS_BUILDING = ['biofuel_boiler', 'electric_heater', 'hp', 'methane_boiler']
HEAT_TECHS_DISTRICT = [
    'chp_biofuel_extraction', 'chp_methane_extraction', 'chp_wte_back_pressure'
]
COOKING = ["electric_hob", "gas_hob"]
VARIABLE_RENEWABLE_TECHS = [
    'open_field_pv',
    'roof_mounted_pv',
    'wind_offshore',
    'wind_onshore',
    "hydro_run_of_river"
]
STORAGE_DISCHARGE_TECHS = [
    'battery_storage',
    'heat_storage_big',
    'heat_storage_small',
    'hydro_storage',
    'hydrogen_storage',
    'ccgt'
]

def calliope_spore_to_friendly_file(
    path_to_result, path_to_annual_demand, path_to_industry_demand, path_to_input_model,
    paths_to_cfs, resolution, year, initial_keywords, name, description, path_to_output,
    scenario, scenario_dim_name
):
    include_hourly = True
    include_ts_in_file = True
    model_dict = {scenario: calliope.read_netcdf(path_to_result)}
    input_model = calliope.read_netcdf(path_to_input_model)
    print("loaded calliope models")

    data_dict = calliope_results_to_friendly(
        model_dict, input_model, path_to_annual_demand, path_to_industry_demand,
        paths_to_cfs, year, scenario, include_hourly, scenario_dim_name
    )

    data_dict = drop_hourly(data_dict)

    keywords = initial_keywords + [
        "calliope", "Sector-coupled Euro-Calliope",
        f"resolution={resolution}H",
        f"weather_year={year}",
        f"model_year=2050"
    ]
    meta = {
        "name": name,
        "description": description,
        "keywords": keywords,
        "licenses": "CC-BY-4.0"
    }
    write_dpkg(
        data_dict, path_to_output, meta, include_timeseries_data=include_ts_in_file
    )


def get_year(path_to_result):
    return int(os.path.basename(path_to_result).split("_")[1])


def calliope_results_to_friendly(
        model_dict, input_model, path_to_annual_demand, path_to_industry_demand,
        paths_to_cfs, year, scenario, include_hourly, scenario_dim_name
    ):
    region_group = None
    data_dict = consolidate_calliope_output.combine_scenarios_to_one_dict(
        model_dict,
        cost_optimal_model=input_model,
        new_dimension_name=scenario_dim_name,
        return_hourly=include_hourly,
        region_group=region_group
    )
    print("loaded base data")
    data_dict_national = consolidate_calliope_output.combine_scenarios_to_one_dict(
        model_dict,
        cost_optimal_model=input_model,
        new_dimension_name=scenario_dim_name,
        return_hourly=include_hourly,
        region_group="countries"
    )
    print("loaded base, nationally aggregated data")
    for variable in ["flow_out", "flow_in", "net_import"]:
        for agg in ["max", "sum"]:
            data_dict[f"{variable}_{agg}_1M"] = data_dict_national[f"{variable}_{agg}_1M"]

    data_dict = add_extra_variables(
        data_dict, path_to_annual_demand, year, scenario, scenario_dim_name,
        path_to_industry_demand, input_model, paths_to_cfs, region_group, model_dict
    )
    print("Added extra variables")

    return data_dict


def add_extra_variables(
    data_dict, path_to_annual_demand, year, scenario, scenario_dim_name,
    path_to_industry_demand, input_model, paths_to_cfs, region_group, model_dict
):
    data_dict["primary_energy_supply"] = add_primary_energy_supply(data_dict)
    print("Added primary energy")

    final_consumption_subnational, final_consumption_industry_subsectors = add_final_energy_consumption(
        data_dict, path_to_annual_demand, year, path_to_industry_demand,
        scenario, scenario_dim_name
    )
    data_dict["final_consumption"] = final_consumption_subnational
    data_dict["industry_subsector_consumption"] = final_consumption_industry_subsectors
    print("Added final consumption")

    data_dict["curtailment"] = add_curtailment(
        paths_to_cfs, data_dict["nameplate_capacity"],
        data_dict["flow_out_sum"], year, scenario_dim_name
    )
    print("Added curtailment")
    data_dict["fuel_overproduction"] = add_fuel_overproduction(
        final_consumption_subnational, data_dict["flow_in_sum"], scenario_dim_name
    )
    print("Added fuel overprod")
    data_dict["average_national_net_imports"] = add_average_national_imports(
        data_dict["net_import"], scenario_dim_name
    )
    print("Added net imports")
    data_dict["ev_charge_correlation"] = add_ev_charge_correlation(
        data_dict["flow_out"], data_dict["flow_in"], scenario_dim_name
    )
    print("Added ev charge correlation")
    data_dict["grid_capacity_expansion"] = get_grid_capacity_expansion(
        input_model.inputs, data_dict["net_transfer_capacity"], region_group
    )
    print("Added grid capacity expansion")
    data_dict["grid_transfer_capacity"] = data_dict.pop("net_transfer_capacity")

    data_dict["paper_metrics"] = get_paper_metrics(
        data_dict, input_model.inputs, scenario, scenario_dim_name
    )
    print("Added paper metrics")

    return data_dict


def add_primary_energy_supply(data_dict):
    def _sum_over_techs(flow, tech_grouping):
        return (
            data_dict[flow]
            .dropna()
            .unstack("techs")
            .groupby(tech_grouping, axis=1).sum(min_count=1)
            .droplevel("carriers", axis=0)
            .rename_axis(columns="carriers")
            .stack()
            .dropna()
            .sum(level=data_dict[flow].index.names.difference(["techs"]))
        )

    flow_out_summed = _sum_over_techs("flow_out_sum", ENERGY_PRODUCERS)

    flow_out_summed = flow_out_summed.append(
        data_dict["net_import_sum"]
        .rename_axis(index={"importing_region": "locs"})
        .rename({"electricity": "Net electricity import"}, level="carriers")
        .droplevel("exporting_region")
        .sum(level=flow_out_summed.index.names)
    )

    return flow_out_summed.sort_index()


def add_final_energy_consumption(
    data_dict, path_to_annual_demand, year, path_to_industry_demand,
    scenario, scenario_dim_name
):

    road_transport_consumption, heat_consumption_building, heat_consumption_district, all_elec_consumption, idx_order = (
        final_energy_consumption_from_data_dict(data_dict)
    )


    air_transport_consumption, shipping_transport_consumption, rail_transport_consumption, industry_total_elec_consumption, industry_all_consumption = final_energy_consumption_from_annual_demand(
        path_to_annual_demand, idx_order, year, scenario, scenario_dim_name
    )

    countries = all_elec_consumption.index.get_level_values("locs").str.split("_", expand=True).levels[0]
    industry_subsector_consumption = final_energy_consumption_from_industry_subsectors(
        path_to_industry_demand, year, countries, scenario, scenario_dim_name
    )

    building_electricity_consumption = (
        all_elec_consumption
        .sub(rail_transport_consumption)
        .sub(industry_total_elec_consumption)
    )
    final_idx_order = [scenario_dim_name, "sector", "subsector", "carriers", "locs", "unit"]
    add_final_levels = add_sector_subsector(idx_order=final_idx_order)
    all_consumption_subnational = pd.concat([
        add_final_levels(road_transport_consumption, sector="Transport", subsector="Road"),
        add_final_levels(heat_consumption_building, sector="Buildings", subsector="Building heat"),
        add_final_levels(heat_consumption_district, sector="Buildings", subsector="District heat"),
        add_final_levels(air_transport_consumption, sector="Transport", subsector="Aviation"),
        add_final_levels(shipping_transport_consumption, sector="Transport", subsector="Shipping"),
        add_final_levels(rail_transport_consumption, sector="Transport", subsector="Rail"),
        add_final_levels(industry_all_consumption, sector="Industry", subsector="Processes and feedstocks"),
        add_final_levels(building_electricity_consumption, sector="Buildings", subsector="Appliances & cooling")
    ]).sort_index()

    industry_subsector_consumption = add_final_levels(
        industry_subsector_consumption, sector="Industry"
    ).sort_index()

    return all_consumption_subnational, industry_subsector_consumption


def final_energy_consumption_from_data_dict(data_dict):
    flow_in = data_dict["flow_in_sum"].dropna()
    flow_out = data_dict["flow_out_sum"].dropna()
    road_transport_techs = flow_out.xs("transport", level="carriers").unstack("techs").columns
    road_transport_consumption = (
        flow_in
        .unstack("techs")
        .loc[:, road_transport_techs]
        .sum(axis=1, min_count=1)
        .dropna()
    )

    heat_consumption_building = (
        flow_in[flow_in.index.get_level_values("techs").isin(HEAT_TECHS_BUILDING)]
        .unstack("techs")
        .sum(axis=1, min_count=1)
    )
    heat_consumption_district = (
        flow_out[
            flow_out.index.get_level_values("techs").isin(HEAT_TECHS_DISTRICT) &
            (flow_out.index.get_level_values("carriers") == "heat")
        ]
        .unstack("techs")
        .sum(axis=1, min_count=1)
    )
    all_elec_consumption = flow_in.xs("demand_elec", level="techs")

    idx_order = [i for i in flow_out.index.names if i != "techs"]
    return (
        road_transport_consumption, heat_consumption_building, heat_consumption_district,
        all_elec_consumption, idx_order
    )


def final_energy_consumption_from_annual_demand(
    path_to_annual_demand, idx_order, year, scenario, scenario_dim_name
):
    annual_demand = util.read_tdf(path_to_annual_demand)

    air_transport_consumption = annual_demand_to_scenario_format(
        annual_demand
        .xs(("industry_demand", "air", year), level=("dataset", "cat_name", "year")),
        idx_order, scenario, scenario_dim_name
    )
    shipping_transport_consumption = annual_demand_to_scenario_format(
        annual_demand
        .xs(("industry_demand", "marine", year), level=("dataset", "cat_name", "year")),
        idx_order, scenario, scenario_dim_name
    )
    rail_transport_consumption = annual_demand_to_scenario_format(
        annual_demand
        .xs(("rail", year), level=("cat_name", "year"))
        .unstack("end_use")
        .loc[:, ["electricity"]]
        .stack()
        .sum(level=["id", "unit", "end_use"], min_count=1),
        idx_order, scenario, scenario_dim_name
    )
    industry_total_elec_consumption = annual_demand_to_scenario_format(
        annual_demand.xs(
            ("industry_demand", "industry", year),
            level=("dataset", "cat_name", "year")
        )
        .unstack("end_use")
        .loc[:, ["electricity"]]
        .stack(),
        idx_order, scenario, scenario_dim_name
    )
    industry_all_consumption = annual_demand_to_scenario_format(
        annual_demand.xs(
            ("industry_demand", "industry", year),
            level=("dataset", "cat_name", "year")
        ),
        idx_order, scenario, scenario_dim_name
    )
    return (
        air_transport_consumption, shipping_transport_consumption,
        rail_transport_consumption, industry_total_elec_consumption,
        industry_all_consumption
    )


def final_energy_consumption_from_industry_subsectors(
    path_to_industry_demand, year, countries, scenario, scenario_dim_name
):
    industry_demand = util.read_tdf(path_to_industry_demand)
    industry_subsector_consumption = (
        industry_demand
        .xs(year, level="year")
        .unstack("country_code")
        .groupby({util.get_alpha2(i): i for i in countries}, axis=1)
        .sum(min_count=1)
        .rename_axis(columns="locs", index={"carrier": "carriers"})
        .stack()
        .drop("space_heat", level="carriers")
        .to_frame(scenario)
        .rename_axis(columns=scenario_dim_name)
        .stack()
        .reorder_levels([scenario_dim_name, "subsector", "carriers", "locs", "unit"])
    )
    return industry_subsector_consumption


def annual_demand_to_scenario_format(annual_demand, idx_order, scenario, scenario_dim_name):
    units = annual_demand.index.get_level_values("unit").unique()
    dfs = []
    for scaled_unit in units:
        scale, unit = units.str.split(" ")[0]
        scale = float(scale)
        dfs.append(
            annual_demand
            .mul(scale)
            .rename({scaled_unit: unit}, level="unit")
            .rename_axis(index={"id": "locs", "end_use": "carriers"})
        )
    _df = pd.concat(dfs)

    return (
        _df
        .to_frame(scenario)
        .rename_axis(columns=scenario_dim_name)
        .stack()
        .sort_index()
        .reorder_levels(idx_order)
    )


def add_sector_subsector(idx_order):
    def _add_sector_subsector(df, **kwargs):
        for dim_name, dim_val in kwargs.items():
            df = df.to_frame(dim_val).rename_axis(columns=dim_name).stack()
        return df.reorder_levels(idx_order).sort_index()
    return _add_sector_subsector


def _cf_from_file(path_to_cf):
    df = pd.read_csv(path_to_cf, index_col=0, parse_dates=True)
    return df.groupby(df.index.year).sum().rename_axis(index="weather_year", columns="locs")

def _cf_tech_name_from_file(path_to_cf):
    return os.path.basename(path_to_cf).replace("capacityfactors-", "").replace("-", "_").replace(".csv", "")


def add_curtailment(paths_to_cfs, nameplate_capacity, flow_out_sum, year, scenario_dim_name):
    _techs = "_pv|wind_|hydro_"
    cf = pd.concat(
        [_cf_from_file(path_to_cf) for path_to_cf in paths_to_cfs],
        keys=[_cf_tech_name_from_file(path_to_cf) for path_to_cf in paths_to_cfs],
        names=["techs"]
    ).rename({
        "rooftop_pv": "roof_mounted_pv",
        "hydro_ror": "hydro_run_of_river",
        "hydro_reservoir_inflow": "hydro_reservoir"
    }, level="techs")
    _cf = cf.xs(year, level="weather_year").stack()
    technical_maximum = nameplate_capacity.filter(regex=_techs).mul(_cf).dropna()

    return technical_maximum.droplevel("unit").sub(flow_out_sum).dropna()


def add_fuel_overproduction(final_consumption_subnational, flow_in_sum, scenario_dim_name):
    actual_industry_demand = (
        final_consumption_subnational
        .groupby(level=["carriers", "locs", scenario_dim_name])
        .sum()
        .unstack(["locs", scenario_dim_name])
        .loc[["diesel", "kerosene", "methane", "methanol"]]
        .sum()
    )
    produced_fuel_for_industry = (
        flow_in_sum
        .unstack("techs")
        .loc[:, [
            "demand_industry_diesel",
            "demand_industry_kerosene",
            "demand_industry_methane",
            "demand_industry_methanol",
        ]]
        .sum(axis=1)
        .groupby(level=["locs", scenario_dim_name])
        .sum()
        .reindex(actual_industry_demand.index)
        .fillna(0)
    )
    return produced_fuel_for_industry.div(actual_industry_demand)


def add_average_national_imports(net_import, scenario_dim_name):
    return net_import.unstack([scenario_dim_name, "unit"]).groupby(
        [_region_to_country, _region_to_country],
        level=["importing_region", "exporting_region"]
    ).sum().where(lambda x: x > 0).mean()


def _region_to_country(region):
    return region.split("_")[0]


def add_ev_charge_correlation(flow_out, flow_in, scenario_dim_name):
    ev_consumption = (
        flow_in
        .filter(regex="_transport_ev")
        .groupby(level=[scenario_dim_name, "timesteps"])
        .sum()
    )

    renewables_production = (
        flow_out
        .unstack("techs")
        .reindex(VARIABLE_RENEWABLE_TECHS, axis=1)
        .sum(axis=1)
        .groupby(level=[scenario_dim_name, "timesteps"])
        .sum()
    )
    to_correlate = pd.concat(
        [ev_consumption, renewables_production], axis=1, keys=["ev", "renewables"]
    )
    return pd.concat(
        [pd.Series(
            {metric: to_correlate.xs(scenario, level=scenario_dim_name).corr(metric).loc["ev", "renewables"]
            for metric in ["spearman", "pearson"]},
            name="correlation_method"
        ) for scenario in to_correlate.index.get_level_values(scenario_dim_name).unique()],
        keys=to_correlate.index.get_level_values(scenario_dim_name).unique(),
        names=[scenario_dim_name, "correlation_method"]
    )


def get_grid_capacity_expansion(inputs, final_capacities, region_group):
    min_capacity = (
        consolidate_calliope_output.clean_series(
            consolidate_calliope_output.map_da(
                inputs.energy_cap_min,
                region_group=region_group,
                keep_demand=False,
                transmission_only=True
            )
        )
        .unstack("locs")
    )
    _from = "exporting_region"
    _to = "importing_region"

    def _rename_remote(series_or_df, level):
        df = series_or_df.reset_index(level)
        if (df[level].str.find(":") > -1).any():
            remote = df[level].str.split(":", expand=True)[1]
        else:
            remote = df[level]
        df[level] = remote
        return df.set_index(level, append=True)

    def _rename_tech(x):
        return "ac_transmission" if x.startswith("ac") else "dc_transmission"

    min_capacity.index = min_capacity.index.str.split(":", expand=True).rename(["techs", _from])
    min_capacity = (
        _rename_remote(min_capacity, level=_from)
        .rename(_rename_tech, level="techs")
        .rename_axis(columns=_to)
        .stack()
        .div(10)  # 0.1 TWh -> TWh
    )
    min_capacity = min_capacity.groupby(level=min_capacity.index.names).sum()

    return final_capacities.sub(min_capacity, fill_value=0)


def get_paper_metrics(data_dict, inputs, scenario, scenario_dim_name):
    metrics = pd.MultiIndex.from_tuples(
        [("curtailment", "percentage", scenario),
        ("electricity_production_gini", "fraction", scenario),
        ("average_national_import", "twh", scenario),
        ("fuel_autarky_gini", "fraction", scenario),
        ("storage_discharge_capacity", "tw", scenario),
        ("transport_electrification", "percentage", scenario),
        ("heat_electrification", "percentage", scenario),
        ("ev_flexibility", "fraction", scenario),
        ("biofuel_utilisation", "percentage", scenario)],
        names=["metric", "unit", scenario_dim_name]
    )
    metric_df = pd.Series(index=metrics)
    metric_df.loc["curtailment"] = (
        data_dict["curtailment"].xs(scenario, level=scenario_dim_name).sum() /
        data_dict["curtailment"].xs(scenario, level=scenario_dim_name)
        .add(data_dict["flow_out_sum"].xs(scenario, level=scenario_dim_name)).dropna().sum()
    )
    metric_df.loc["electricity_production_gini"] = get_gini(
        data_dict["flow_out_sum"]
        .xs("electricity", level="carriers")
        .xs(scenario, level=scenario_dim_name)
        .unstack("techs")
        .reindex(list(ENERGY_PRODUCERS.keys()) + HEAT_TECHS_DISTRICT, axis=1)
        .sum(axis=1)
        .groupby(level="locs")
        .sum()
    )
    metric_df.loc["average_national_import"] = (
        data_dict["average_national_net_imports"]
        .xs(scenario, level=scenario_dim_name)
        .item()
    )
    metric_df.loc["fuel_autarky_gini"] = get_gini(data_dict["fuel_overproduction"])
    metric_df.loc["storage_discharge_capacity"] = (
        data_dict["nameplate_capacity"]
        .xs(scenario, level=scenario_dim_name)
        .unstack("techs")
        .reindex(STORAGE_DISCHARGE_TECHS, axis=1)
        .sum(axis=1)
        .sum()
    )
    metric_df.loc["transport_electrification"] = (
        data_dict["flow_out_sum"]
        .xs(scenario, level=scenario_dim_name)
        .unstack("techs")
        .reindex(["light_transport_ev", "heavy_transport_ev"], axis=1)
        .sum(axis=1, min_count=1)
        .sum()
        /
        data_dict["flow_out_sum"]
        .xs("transport", level="carriers")
        .xs(scenario, level=scenario_dim_name)
        .sum()
    )
    metric_df.loc["heat_electrification"] = (
        data_dict["flow_out_sum"]
        .xs(scenario, level=scenario_dim_name)
        .unstack("techs")
        .reindex(["hp", "electric_heater", "electric_hob"], axis=1)
        .sum(axis=1)
        .sum()
        / (
            data_dict["flow_out_sum"]
            .xs(scenario, level=scenario_dim_name)
            .unstack("techs")
            .reindex(HEAT_TECHS_BUILDING + HEAT_TECHS_DISTRICT + COOKING , axis=1)
            .sum(axis=1)
            .unstack("carriers")
            .loc[:, ["heat", "cooking"]]
            .sum(axis=1)
            .sum()
        )
    )
    metric_df.loc["ev_flexibility"] = (
        data_dict["ev_charge_correlation"]
        .xs((scenario, "pearson"), level=(scenario_dim_name, "correlation_method"))
        .item()
    )
    if "biofuel_supply" in data_dict["flow_out_sum"].xs(scenario, level=scenario_dim_name).index.get_level_values("techs"):
        metric_df.loc["biofuel_utilisation"] = (
            data_dict["flow_out_sum"]
            .xs(scenario, level=scenario_dim_name)
            .xs("biofuel_supply", level="techs")
            .sum()
            / (
                inputs
                .group_carrier_prod_max
                .to_series()
                .filter(regex="biofuel")
                .div(10)  # 0.1 TWh -> TW
                .sum()
            )
        )
    else:
        metric_df.loc["biofuel_utilisation"] = 0
    metric_df.loc[metric_df.index.get_level_values("unit") == "percentage"] *= 100
    return metric_df


def get_gini(metric):
    """
    Get the gini index for a particular metric.
    This is used to give an indication of the spatial 'equity' of a metric.
    """

    results = []
    vals = metric.values
    for i in vals:
        for j in vals:
            results.append(abs(i - j))
    return sum(results) / (2 * len(vals)**2 * vals.mean())


def drop_hourly(data_dict):
    for hourly_metric in ["flow_in", "flow_out", "net_import", "storage_level"]:
        del data_dict[hourly_metric]

    return data_dict


if __name__ == "__main__":
    calliope_spore_to_friendly_file(
        path_to_result=snakemake.input.model_result,
        path_to_annual_demand=snakemake.input.annual_demand,
        path_to_industry_demand=snakemake.input.industry_demand,
        path_to_input_model=snakemake.input.input_model,
        paths_to_cfs=snakemake.input.cfs,
        initial_keywords=snakemake.params.initial_keywords,
        name=snakemake.params.name,
        description=snakemake.params.description,
        resolution=snakemake.wildcards.get("model_resolution", 2),
        year=int(snakemake.params.year),
        scenario=snakemake.params.scenario,
        scenario_dim_name=snakemake.params.scenario_dim_name,
        path_to_output=snakemake.output[0]
    )
