import sys

import pandas as pd
import geopandas as gpd
import jinja2

import util

sys.path.append('euro-calliope/src/')
import filters


# TODO get charge/disharge efficiency
TEMPLATE = """
overrides:
    gas_storage:
        techs.methane_storage:
            essentials:
                name: Underground methane storage
                parent: storage
                carrier: methane
            constraints:
                energy_eff: 1

        locations:
            {% for idx in gas_storage.index %}
            {% if gas_storage.loc[idx].sum() > 0 %}
            {{ idx }}:
                techs.methane_storage.constraints:
                    energy_cap_equals: {{ gas_storage.loc[idx, 'energy_cap_gas_storage_mw'] * scaling_factors.power }}  # {{ (1 / scaling_factors.power) | unit("MW") }}
                    storage_cap_equals: {{ gas_storage.loc[idx, 'storage_cap_gas_storage_mwh'] * scaling_factors.power }}  # {{ (1 / scaling_factors.power) | unit("MWh") }}
            {% endif %}
            {% endfor %}
"""


def regionalise_gas_storage(
    path_to_storage_data, path_to_units, path_to_table_results, path_to_yaml_results,
    scaling_factors
):
    """
    Take underground natural gas storage data from Gas Infrastructure Europe
    and assign it to model regions.

    Data comes from an Excel spreadsheet and requires a bit of massaging to
    handle some non-machine-readable elements.

    Currently assigning storage info evenly to all regions on a subnational level.
    Would be better to use exact locatonal data, as shown here:
    https://www.gie.eu/download/maps/2018/GIE_STOR_2018_A0_1189x841_FULL_FINAL.pdf
    """

    all_sites = pd.read_excel(
        path_to_storage_data, sheet_name="Storage DB ", skiprows=1, index_col=0
    )

    daily_energy_cap = (
        all_sites
        .loc[(all_sites.index == 'Serbia') | (all_sites["in EU28 SUM"] == "y"),
             ["Withdrawal\ntechnical\nGWh/day", "Injection\ntechnical\nGWh/day"]]
        .where(lambda x: x > 0)
        .mean(axis=1)
    )
    hourly_energy_cap = (daily_energy_cap / 24).groupby(level=0).sum()

    gas_storage_cap = pd.read_excel(
        path_to_storage_data, sheet_name="WGV by type TWh", skiprows=3, index_col=0
    )

    gas_storage_data = pd.concat(
        [hourly_energy_cap * 1e3,  # Scale GW to MW
         gas_storage_cap.reindex(hourly_energy_cap.index).Total * 1e6],  # scale TWh to MWh
        keys=['energy_cap_gas_storage_mw', 'storage_cap_gas_storage_mwh'],
        axis=1
    )
    gas_storage_data.index = gas_storage_data.index.map(util.get_alpha3).rename('country_code')

    units = gpd.read_file(path_to_units).set_index(["id", "country_code"])

    # TODO: split gas to subcountry levels in a way that is based on actual postion of
    # storage facilities. Currently being split evenly across all regions.
    gas_storage_data_country = gas_storage_data.reindex(units.index, level="country_code")
    gas_storage_data_per_region = (
        gas_storage_data_country
        .div(gas_storage_data_country.count(level='country_code'))
        .droplevel('country_code')
        .fillna(0)
    )

    gas_storage_data_per_region.to_csv(path_to_table_results)


    env = jinja2.Environment(lstrip_blocks=True, trim_blocks=True)
    env.filters["unit"] = filters.unit

    gas_storage_yaml = env.from_string(TEMPLATE).render(
        gas_storage=gas_storage_data_per_region,
        scaling_factors=scaling_factors,
    )

    with open(path_to_yaml_results, "w") as result_file:
        result_file.write(gas_storage_yaml)


if __name__ == "__main__":
    regionalise_gas_storage(
        path_to_storage_data=snakemake.input.gas_storage_data,
        path_to_units=snakemake.input.units,
        scaling_factors=snakemake.params.scaling_factors,
        path_to_table_results=snakemake.output.table,
        path_to_yaml_results=snakemake.output.yaml
    )
