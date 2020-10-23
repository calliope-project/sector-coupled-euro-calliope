import sys

import jinja2
import pandas as pd
import pycountry

import util

sys.path.append('euro-calliope/src/')
import filters

TEMPLATE = """
overrides:
    annual_fuel_demand:
        group_constraints:
            {% for idx in annual_demand.index %}
            {% for demand, demand_key in demand_keys.items() %}
            {{ demand }}_{{ idx }}:
                techs: [{{ demand }}]
                locs: [{{ idx }}]
                carrier_con_equals:
                    {{ demand_key[1] }}: {{ annual_demand.loc[idx, demand_key] }}  # {{ demand_key[2] }}
            {% endfor %}

            {% endfor %}
    industry_techs:
        locations:
        {% for idx in annual_demand.index %}
            {{ idx }}.techs:
                biofuel_to_liquids:
                hydrogen_to_liquids:
                biofuel_to_diesel:
                biofuel_to_methanol:
                hydrogen_to_methanol:
                biofuel_to_methane:
                hydrogen_to_methane:
                electrolysis:
                dac:
                demand_air_kerosene:
                demand_marine_diesel:
                demand_industry_hydrogen:
                demand_industry_methanol:
                demand_industry_co2:
                demand_industry_methane:
                demand_industry_diesel:
        {% endfor %}
    new_hydrogen_storage:
        techs:
            hydrogen: # based on [Danish energy agency, energy storage, 151a Hydrogen Storage - Tanks, 2050]
                essentials:
                    name: Hydrogen power storage
                    parent: storage
                    carrier: hydrogen
                constraints:
                    energy_cap_max: inf
                    storage_cap_max: inf
                    energy_cap_per_storage_cap: 0.0048
                    energy_eff: 0.81  # 0.90 round-trip
                    lifetime: 30
                costs:
                    monetary:
                        storage_cap: {{ 0.021e6 * scaling_factors.specific_costs }}  # {{ (1 / scaling_factors.specific_costs) | unit("EUR2015/MWh_hydrogen") }}
                        om_annual: {{ 400 * scaling_factors.specific_costs }} # {{ (1 / scaling_factors.specific_costs) | unit("EUR2015/MW_hydrogen/year") }}
    new_biofuel_supply:
        techs:
            biofuel_supply:
                essentials:
                    name: Biofuel
                    parent: supply
                    carrier: biofuel
                # constraints are applied on a per-location basis
                costs.monetary:
                    energy_cap: 0  # negate the costs given in the base euro-calliope definition (anaerobic digection gas turbine)
                    om_prod: 0  # negate the costs given in the base euro-calliope definition (anaerobic digection gas turbine)
                    om_annual: 0  # negate the costs given in the base euro-calliope definition (anaerobic digection gas turbine)
                    om_con: {{ biofuel_fuel_cost * scaling_factors.specific_costs }} # {{ (1 / scaling_factors.specific_costs) | unit("EUR/MWh") }}
scenarios:
    industry_fuel: [annual_fuel_demand, industry_techs, new_hydrogen_storage, new_biofuel_supply, biofuel_maximum]
"""


def fill_constraint_template(path_to_annual_demand, path_to_result, model_year, scaling_factors, path_to_biofuel_costs):
    """Generate a file that represents links in Calliope."""
    annual_demand = util.read_tdf(path_to_annual_demand)
    annual_demand = annual_demand.xs(('industry_demand', model_year), level=('dataset', 'year'))
    keys = (
        annual_demand
        .droplevel(['id'])
        .index.drop_duplicates()
        .drop(['rail', 'road'], level='cat_name')
        .drop(['bau_electricity', 'electricity', 'space_heat'], level='end_use')
        .reorder_levels(['cat_name', 'end_use', 'unit'])
    )
    demand_keys = {'demand_{}_{}'.format(*i[:-1]): i for i in keys}
    scaling_factors["specific_costs"] = scaling_factors["monetary"] / scaling_factors["power"]
    with open(path_to_biofuel_costs, "r") as f_biofuel_costs:
        biofuel_fuel_cost = float(f_biofuel_costs.readline())

    env = jinja2.Environment(lstrip_blocks=True, trim_blocks=True)
    env.filters["unit"] = filters.unit

    links = env.from_string(TEMPLATE).render(
        annual_demand=annual_demand.unstack(keys.names)[keys].fillna(0),
        demand_keys=demand_keys,
        scaling_factors=scaling_factors,
        biofuel_fuel_cost=biofuel_fuel_cost
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(links)


if __name__ == "__main__":
    fill_constraint_template(
        path_to_annual_demand=snakemake.input.annual_demand,
        path_to_result=snakemake.output[0],
        model_year=snakemake.params.model_year,
        scaling_factors=snakemake.params.scaling_factors,
        path_to_biofuel_costs=snakemake.input.biofuel_cost
    )
