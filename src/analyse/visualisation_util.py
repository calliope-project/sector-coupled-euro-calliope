import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import calliope
from calliope.core.util.dataset import split_loc_techs
idx = pd.IndexSlice

CONFIG = calliope.AttrDict.from_yaml("../config/default.yaml")


class VisUtil():

    VARIABLE_RENEWABLE_TECHS = [
        'open_field_pv',
        'roof_mounted_pv',
        'wind_offshore',
        'wind_onshore_competing',
        'wind_onshore_monopoly',
    ]
    FIRM_SUPPLY_TECHS = [
        "ccgt", "chp_biofuel_extraction",
        "chp_methane_extraction", "chp_wte_back_pressure",
        "nuclear", "hydro_reservoir", "hydro_run_of_river"
    ]
    ELECTRICITY_SUPPLY_TECHS = [
        "ccgt", "chp_biofuel_extraction",
        "chp_methane_extraction", "chp_wte_back_pressure",
        "nuclear", "open_field_pv", "roof_mounted_pv", "wind_offshore",
        "wind_onshore_competing", "wind_onshore_monopoly", "hydro_reservoir",
        "hydro_run_of_river"
    ]
    RENEWABLES_MAPPING = {
        'open_field_pv': 'pv',
        'roof_mounted_pv': 'pv',
        'wind_offshore': 'offshore wind',
        'wind_onshore_competing': 'onshore wind',
        'wind_onshore_monopoly': 'onshore wind',
        'biofuel_supply': 'biofuel_supply'
    }
    FIRM_SUPPLY_MAPPING = {
        "ccgt": "ccgt",
        "chp_biofuel_extraction": "chp",
        "chp_methane_extraction": "chp",
        "chp_wte_back_pressure": "chp",
        "nuclear": "nuclear",
        "hydro_reservoir": "hydro",
        "hydro_run_of_river": "hydro",
    }
    STORAGE_MAPPING = {
        "hp_heat_storage_small": "heat_small",
        "electric_heater_heat_storage_small": "heat_small",
        "biofuel_heat_storage_small": "heat_small",
        "methane_heat_storage_small": "heat_small",
        "chp_wte_back_pressure_heat_storage_big": "heat_big",
        "chp_biofuel_extraction_heat_storage_big": "heat_big",
        "chp_methane_extraction_heat_storage_big": "heat_big",
        "methane_storage": "gas_storage" ,
        "hydrogen_storage": "gas_storage" ,
        "pumped_hydro": "electricity_storage",
        "battery": "electricity_storage",
    }
    TRANSMISSION_MAPPING = {
        'ac_ohl_mountain_transmission': 'ac_transmission',
        'ac_ohl_transmission': 'ac_transmission',
        'dc_ohl_transmission': 'dc_transmission',
        'dc_subsea_transmission': 'dc_transmission',
        'dc_underground_transmission': 'dc_transmission',
    }
    HEAT_TECH_MAPPING = {
        "electric_heater": "direct",
        "hp": "direct",
        "chp_biofuel_extraction": "district",
        "chp_wte_back_pressure": "district",
        "biofuel_boiler": "direct",
        "methane_boiler": "direct",
        "chp_methane_extraction": "district",
    }
    HEAT_TECH_CARRIER_MAPPING = {
        "electric_heater": "electric",
        "hp": "electric",
        "electric_hob": "electric",
        "chp_biofuel_extraction": "biofuel",
        "chp_wte_back_pressure": "waste",
        "biofuel_boiler": "biofuel",
        "methane_boiler": "methane",
        "chp_methane_extraction": "methane",
        "gas_hob": "methane"
    }
    TRANSPORT_MAPPING = {
        "heavy_transport_ice": "ice",
        "light_transport_ice": "ice",
        "heavy_transport_ev": "ev",
        "light_transport_ev": "ev",
    }
    SYNTHETIC_FUEL_MAPPING = {
        "biofuel_to_liquids": "biofuel",
        "biofuel_to_diesel": "biofuel",
        "biofuel_to_gas_and_liquids": "biofuel",
        "biofuel_to_methanol": "biofuel",
        "biofuel_to_methane": "biofuel",
        "hydrogen_to_liquids": "hydrogen",
        "hydrogen_to_methanol": "hydrogen",
        "hydrogen_to_methane": "hydrogen",
    }


    def __init__(self, scenario_name, model, inputs=None):
        self.scenario = scenario_name
        self.model_data = model._model_data
        if inputs is None:
            self.inputs = self.model_data.filter_by_attrs(is_result=0)
        else:
            self.inputs = inputs
        self.set_power_density(CONFIG)

        self.carrier_prod = self.clean_series(self.model_data.carrier_prod.sum("timesteps", min_count=1))
        self.energy_cap = self.clean_series(self.model_data.energy_cap)
        self.storage_cap = self.clean_series(self.model_data.storage_cap)
        self.energy_cap_max = self.clean_series(self.inputs.energy_cap_max)
        self.group_carrier_prod_max = self.clean_series(self.inputs.group_carrier_prod_max)
        self.available_area = self.inputs.available_area.to_series()

        self.metrics = {}
        self.set_all_metrics()

    def set_power_density(self, config):
        power_density = config["parameters"]['maximum-installable-power-density'].copy()  # MW/km^2
        power_density['wind_onshore_monopoly'] = power_density['onshore-wind']
        power_density['wind_onshore_competing'] = power_density.pop('onshore-wind')
        power_density['wind_offshore'] = power_density.pop('offshore-wind')
        power_density['roof_mounted_pv'] = power_density.pop('pv-on-tilted-roofs')
        power_density['open_field_pv'] = power_density.pop('pv-on-flat-areas')

        self.power_density = pd.Series({
            k: v * config["scaling-factors"]["power"] / config["scaling-factors"]["area"]
            for k, v in power_density.items()
        })

    @staticmethod
    def clean_series(da):

        initial_series = split_loc_techs(da, return_as="Series").where(lambda x: ~np.isinf(x)).dropna()
        initial_levels = initial_series.index.names
        if "techs" in initial_levels:
            non_tech_levels = [i for i in initial_series.index.names if i != "techs"]
            tech_focussed_series = initial_series.unstack(non_tech_levels)
            tech_focussed_series.index = tech_focussed_series.index.str.split(':', expand=True)
            return (
                tech_focussed_series
                .groupby(level=0).sum(min_count=1)
                .rename_axis(index="techs")
                .unstack()
                .reorder_levels(initial_levels)
                .dropna()
            )

        else:
            return initial_series

    @staticmethod
    def add_tech_level_to_series(series, tech_name):
        return series.to_frame(tech_name).rename_axis(columns="techs").stack()

    def sum_then_groupby(self, series, techs, keep_locs=False):
        levels_to_keep, _idx = self.get_levels_and_idx(keep_locs, techs)
        return (
            series
            .rename(techs, level="techs")
            .groupby(level=levels_to_keep).sum()
            .loc[_idx]
        )

    def just_sum(self, series, techs, keep_locs=False):
        levels_to_keep, _idx = self.get_levels_and_idx(keep_locs, techs)
        return series.sum(level=levels_to_keep).loc[_idx]

    @staticmethod
    def get_levels_and_idx(keep_locs, techs):
        try:
            _techs = list(set(techs.values()))
        except AttributeError:
            _techs = techs
        if keep_locs:
            levels_to_keep = ["locs", "techs"]
            _idx = idx[:, _techs]
        else:
            levels_to_keep = "techs"
            _idx = _techs
        return levels_to_keep, _idx

    def clean_transmission(self):
        flows = []
        for flow in ["prod", "con"]:
            _flow = (
                split_loc_techs(
                    self.model_data[f"carrier_{flow}"].sum("timesteps", min_count=1), return_as="Series"
                )
                .where(lambda x: ~np.isinf(x))
                .dropna()
                .filter(regex="transmission")
                .xs("electricity")
                .unstack('locs')
            )
            remote_loc = "loc_to" if flow == "con" else "loc_from"
            loc = "loc_from" if flow == "con" else "loc_to"
            _flow.index = _flow.index.str.split(":", expand=True).rename(["techs", remote_loc])
            _flow = _flow.rename_axis(columns=loc).stack().sum(level=["loc_from", "loc_to"])
            flows.append(_flow.abs())
        transmission = pd.concat(flows, axis=1).mean(axis=1)
        transmission = pd.concat(
            [transmission, transmission.rename_axis(index=["loc_to", "loc_from"]).reorder_levels([1, 0])],
            axis=1, keys=["import", "export"]
        )
        # negative == net export from loc_from to loc_to
        # positive == net import to loc_from from loc_to
        return transmission["import"] - transmission["export"]

    def set_resource_availability(self):
        energy_cap_max = self.energy_cap_max.loc[idx[:, self.power_density.index]]
        biofuel_availability = self.group_carrier_prod_max
        biofuel_availability.index = biofuel_availability.index.str.split("_", 1, expand=True).rename(["techs", "locs"]).reorder_levels(["locs", "techs"])
        biofuel_availability = biofuel_availability.rename(index={"biofuel": "biofuel_supply"})
        available_area = pd.concat([
            self.add_tech_level_to_series(self.available_area, "open_field_pv"),
            self.add_tech_level_to_series(self.available_area, "wind_onshore_competing")
        ]).mul(self.power_density, level="techs")

        self.metrics["resource_availability"] = pd.concat([
            energy_cap_max,
            biofuel_availability,
            available_area
        ])

    def set_resource_use(self):
        if "resource_availability" not in self.metrics.keys():
            self.set_resource_availability()
        energy_cap = self.energy_cap.reorder_levels(self.metrics["resource_availability"].index.names).reindex(self.metrics["resource_availability"].index).fillna(0)
        energy_cap.update(self.carrier_prod.xs("biofuel", level="carriers"))

        self.metrics["resource_use"] = energy_cap

    def set_resource_use_share(self):
        if "resource_availability" not in self.metrics.keys():
            self.set_resource_availability()
        if "resource_use" not in self.metrics.keys():
            self.set_resource_use()

        resource_use_grouped = self.sum_then_groupby(self.metrics["resource_use"], self.RENEWABLES_MAPPING, keep_locs=True)
        availability_grouped = self.sum_then_groupby(self.metrics["resource_availability"], self.RENEWABLES_MAPPING, keep_locs=True)
        self.metrics["resource_use_share"] = resource_use_grouped.div(availability_grouped)
        self.metrics["resource_use_share_total"] = resource_use_grouped.sum(level="techs").div(availability_grouped.sum(level="techs"))

    def set_pv_vs_wind(self):
        self.metrics["pv_vs_wind_cap"] = self.sum_then_groupby(self.energy_cap, self.RENEWABLES_MAPPING, keep_locs=True)
        self.metrics["pv_vs_wind_prod"] = self.sum_then_groupby(self.carrier_prod, self.RENEWABLES_MAPPING, keep_locs=True)

    def set_clean_and_firm(self):
        self.metrics["clean_and_firm_cap"] = self.sum_then_groupby(self.energy_cap, self.FIRM_SUPPLY_MAPPING, keep_locs=True)
        self.metrics["clean_and_firm_prod"] = self.sum_then_groupby(self.carrier_prod, self.FIRM_SUPPLY_MAPPING, keep_locs=True)

    def set_electricity_total(self):
        self.metrics["electricity_total_cap"] = self.just_sum(self.energy_cap, self.ELECTRICITY_SUPPLY_TECHS)
        self.metrics["electricity_total_prod"] = self.just_sum(self.carrier_prod.xs("electricity", level="carriers"), self.ELECTRICITY_SUPPLY_TECHS)

    def set_storage_requirements(self):
        self.metrics["storage_energy_cap"] = self.sum_then_groupby(self.energy_cap, self.STORAGE_MAPPING)
        self.metrics["storage_storage_cap"] = self.sum_then_groupby(self.storage_cap, self.STORAGE_MAPPING)
        self.metrics["storage_prod"] = self.sum_then_groupby(self.carrier_prod, self.STORAGE_MAPPING)

    def set_ac_vs_dc(self):
        self.metrics["ac_vs_dc_energy"] = self.sum_then_groupby(self.energy_cap.div(2), self.TRANSMISSION_MAPPING)
        self.metrics["ac_vs_dc_prod"] = self.sum_then_groupby(self.carrier_prod, self.TRANSMISSION_MAPPING)

    def set_ac_vs_dc_net_transmission(self):
        self.metrics["transmission_flows"] = self.clean_transmission()

    def set_system_cost(self):
        self.metrics["system_cost"] = self.model_data.cost.loc[{"costs": "monetary"}].sum().item()

    def set_direct_vs_district(self):
        self.metrics["direct_vs_district_energy"] = self.sum_then_groupby(self.energy_cap, self.HEAT_TECH_MAPPING)
        self.metrics["direct_vs_district_prod"] = self.sum_then_groupby(self.carrier_prod, self.HEAT_TECH_MAPPING)

    def set_dependence_on_electrolysis(self):
        self.metrics["electrolysis_energy"] = self.just_sum(self.energy_cap, ["electrolysis"], keep_locs=True)
        self.metrics["electrolysis_prod"] = self.just_sum(self.carrier_prod, ["electrolysis"], keep_locs=True)

    def set_heat_fuel_source(self):
        self.metrics["heat_prod"] = self.sum_then_groupby(self.carrier_prod, self.HEAT_TECH_CARRIER_MAPPING)

    def set_transport_fuel_source(self):
        self.metrics["transport_prod"] = self.sum_then_groupby(self.carrier_prod, self.TRANSPORT_MAPPING)

    def set_synthetic_fuel_source(self):
        self.metrics["synfuel_prod"] = self.sum_then_groupby(self.carrier_prod, self.SYNTHETIC_FUEL_MAPPING)


    def set_all_metrics(self):
        self.set_resource_use_share()
        self.set_electricity_total()
        self.set_pv_vs_wind()
        self.set_clean_and_firm()
        self.set_storage_requirements()
        self.set_ac_vs_dc()
        self.set_system_cost()
        self.set_direct_vs_district()
        self.set_dependence_on_electrolysis()
        self.set_ac_vs_dc_net_transmission()
        self.set_heat_fuel_source()
        self.set_transport_fuel_source()
        self.set_synthetic_fuel_source()


def get_grouped_metrics(scenario_utils, scenario_name):
    grouped_metrics = {}
    for metric in scenario_utils[0].metrics.keys():
        if not isinstance(scenario_utils[0].metrics[metric], pd.Series):
            grouped_metrics[metric] = pd.Series(
                data=[i.metrics[metric] for i in scenario_utils],
                index=pd.Index([i.scenario for i in scenario_utils], name=scenario_name),
            )
        else:
            grouped_metrics[metric] = pd.concat(
                [i.metrics[metric] for i in scenario_utils],
                keys=[i.scenario for i in scenario_utils],
                names=scenario_name,
                axis=1
            )
        grouped_metrics[metric].name = metric
    return grouped_metrics


def df_to_plot(df, name=None):
    if name is None:
        name = df.name
    return df.stack().to_frame(name).reset_index()


def plot_renewables_contribution_to_annual_electricity_prod(grouped_metrics, spores):
    with sns.plotting_context("paper", font_scale=1.2):
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        ylabel = 'Contribution to annual electricity production %'
        df = df_to_plot(
            grouped_metrics["pv_vs_wind_prod"].sum(level="techs").drop("biofuel_supply").div(grouped_metrics["electricity_total_prod"].sum()).mul(100),
            name=ylabel
        )
        ax = sns.lineplot(
            data=df[df.spores.isin(spores)],
            x='techs', y=ylabel, hue="spores",
            marker="o", lw="1", zorder=10, markeredgecolor=None,
            palette=["#03a9fc80" if i is not 0 else "#fc032480" for i in spores],
            legend=False, ax=ax
        )
        sns.violinplot(
            data=df,
            x='techs', y=ylabel, inner="box", color="grey",
            cut=0, ax=ax
        )
        for i in ax.collections:
            i.set_alpha(0.5)

        sns.despine(ax=ax)
    return fig, ax


def plot_renewables_production_regionalisation(grouped_metrics, spores):
    with sns.plotting_context("paper", font_scale=1.2):
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        ylabel = 'Relative difference in maximum and minimum\nregional contribution to per-technology electricity production'
        df = df_to_plot(
            grouped_metrics["pv_vs_wind_prod"].div(grouped_metrics["pv_vs_wind_prod"].sum(level="techs")).unstack('techs').agg(lambda x: x.max() - x.min()).unstack('spores').drop('biofuel_supply'),
            name=ylabel
        )
        ax = sns.lineplot(
            data=df[df.spores.isin(spores)],
            x='techs', y=ylabel, hue="spores",
            marker="o", lw="1", zorder=10, markeredgecolor=None,
            palette=["#03a9fc80" if i is not 0 else "#fc032480" for i in spores],
            legend=False, ax=ax
        )
        sns.violinplot(
            data=df,
            x='techs', y=ylabel, inner="box", color="grey",
            cut=0, ax=ax
        )
        for i in ax.collections:
            i.set_alpha(0.5)

        sns.despine(ax=ax)
    return fig, ax

