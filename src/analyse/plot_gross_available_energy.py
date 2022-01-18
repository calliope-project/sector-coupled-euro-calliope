import sys
import os

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from friendly_data.converters import to_df
from frictionless.package import Package

sys.path.append(os.getcwd())
from src.construct import util

ENERGY_BALANCE_GROUPS = {
    'E7000': 'Electricity',
    'H8000': 'Direct heat',
    'C0000X0350-0370': "Other fossils",
    "C0350-0370": "Other fossils",
    "P1000": "Other fossils",
    'G3000': "Natural gas",
    'W6100_6220': "Waste",
    'N900H': 'Nuclear',
    '^O4000.*$': "Oil",
    'S2000': "Oil",
    "^RA1.*$|^RA2.*$|^RA3.*$|^RA4.*$|^RA5.*$": 'Renewables',
    "^R5.*$|W6210": 'Biofuels',
}

ENERGY_PRODUCERS = {
    "waste_supply": "Waste",
    "biofuel_supply": "Biofuels",
    "hydro_reservoir": "Renewables",
    "hydro_run_of_river": "Renewables",
    "nuclear": "Nuclear",
    "open_field_pv": "Renewables",
    "roof_mounted_pv": "Renewables",
    "wind_offshore": "Renewables",
    "wind_onshore": "Renewables",
}

COLORS = {
    "Oil": "#5d5d5d",
    "Natural gas": "#b9b9b9",
    "Other fossils": "#181818",
    "Nuclear": "#cc0000",
    "Biofuels": "#8fce00",
    "Renewables": "#2986cc",
    "Waste": "#ce7e00",
    "Direct heat": "#f6b26b",
    "Electricity": "#2986cc"
}

NUCLEAR_HEAT_MUTIPLIER = 3  # Eurostat uses a multiplier of 3 to go from nuclear power to nuclear heat


def plot_energy_bars(
    path_to_energy_balances, path_to_spores, countries, model_year, path_to_output
):
    country_subselection, input_energy = get_input_energy(
        path_to_energy_balances, countries, model_year
    )
    smallest_output_energy, biggest_spore_energy = get_output_energy(
        path_to_spores, country_subselection
    )
    all_data = pd.concat(
        [input_energy, smallest_output_energy, biggest_spore_energy],
        keys=[f"{model_year:.0f} actual", "Lowest energy\nSPORE", "Highest energy\nSPORE"],
        axis=1
    ).sort_values(f"{model_year:.0f} actual", ascending=False)

    all_data_to_plot = (
        all_data
        .div(1000)
        .where(all_data.div(all_data.sum()) > 0.005)
        .dropna(how="all")
        .T
    )
    with sns.plotting_context("paper", font_scale=2):
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))

        all_data_to_plot.plot.bar(
            ax=ax, stacked=True,
            color=[COLORS[i] for i in all_data_to_plot.columns]
        )
        for tick in ax.get_xticklabels():
            tick.set_rotation(0)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(
            handles[::-1], labels[::-1],
            bbox_to_anchor=(1, 0.5),
            frameon=False,
            loc="center left"
        )

        sns.despine(ax=ax)
        ax.set_ylabel("1000 TWh")

        if path_to_output.endswith(".png"):
            kwargs = {"dpi": 300}
        else:
            kwargs = {}

        fig.savefig(path_to_output, bbox_inches="tight", **kwargs)


def get_input_energy(path_to_energy_balances, countries, model_year):
    annual_energy_balances = util.read_tdf(path_to_energy_balances)

    countries_alpha2 = [util.get_alpha2(i) for i in countries]

    gross_avail_energy = (
        annual_energy_balances
        .xs(("GAE", "TJ", model_year), level=("cat_code", "unit", "year"))
        .unstack()
        .reindex(countries_alpha2, axis=1)
    )
    missing_countries = gross_avail_energy.columns[gross_avail_energy.isnull().all()]
    if not missing_countries.empty:
        print(f"Missing {missing_countries} from input gross annual energy balances")

    to_group = (
        gross_avail_energy
        .index
        .to_frame()
        .squeeze()
        .replace(ENERGY_BALANCE_GROUPS, regex=True)
    )
    to_group = to_group[to_group.index != to_group.values]

    grouped_gross_avail_energy = (
        gross_avail_energy
        .groupby(to_group).sum()
        .sum(axis=1)
        .sort_values(ascending=False)
        .apply(util.tj_to_twh)
    )
    grouped_gross_avail_energy.loc["Nuclear"] /= NUCLEAR_HEAT_MUTIPLIER
    country_subselection = gross_avail_energy.dropna(how="all", axis=1).columns
    return country_subselection, grouped_gross_avail_energy


def get_output_energy(path_to_spores_dpkg, country_subselection):
    spores_datapackage = Package(
        os.path.join(path_to_spores_dpkg, "spores", "datapackage.json")
    )
    spores_data_dicts = {
        resource["name"]: to_df(resource).squeeze().unstack("scenario")
        for resource in spores_datapackage["resources"]
    }
    countries = "|".join([util.get_alpha3(i) for i in country_subselection])
    flow_out = spores_data_dicts["flow_out_sum"]
    flow_out_summed = (
        flow_out[flow_out.index.get_level_values("region").str.contains(countries)]
        .sum(level=["technology", "carriers"])
        .groupby(ENERGY_PRODUCERS, level="technology").sum()
    )
    smallest_spore = flow_out_summed.sum().idxmin()
    biggest_spore = flow_out_summed.sum().idxmax()
    return flow_out_summed[smallest_spore], flow_out_summed[biggest_spore]


if __name__ == "__main__":
    plot_energy_bars(
        path_to_energy_balances=snakemake.input.energy_balances,
        path_to_spores=snakemake.input.spores,
        countries=snakemake.params.countries,
        model_year=snakemake.params.model_year,
        path_to_output=snakemake.output[0]
    )
