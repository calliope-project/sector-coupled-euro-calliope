import sys

import jinja2
import pandas as pd
import pycountry

import util

sys.path.append('euro-calliope/src/')
import filters


TEMPLATE = """
overrides:
    annual_fuel_demand_isolated:
        group_constraints:
            {% for idx in annual_demand.index %}
            {% for carrier in carriers %}
            {% if annual_demand.loc[idx, carrier] > 1e-6 %}
            {{ idx }}_{{ carrier }}:
                techs: [demand_industry_{{ carrier }}]
                locs: [{{ idx }}]
                carrier_con_min:
                    {{ carrier }}: {{ -1 * annual_demand.loc[idx, carrier] * timedelta }}
            {% endif %}
            {% endfor %}
            {% endfor %}

    annual_fuel_demand_shared:
        group_constraints:
            {% for carrier in carriers %}
            {% if annual_demand[carrier].sum() > 1e-6 %}
            all_{{ carrier }}:
                techs: [demand_industry_{{ carrier }}]
                carrier_con_min:
                    {{ carrier }}: {{ -1 * annual_demand[carrier].sum() * timedelta }}
            {% endif %}
            {% endfor %}

    industry_techs:
        techs:
            {% for carrier in carriers %}
            demand_industry_{{ carrier }}.exists: true
            {% endfor %}

        locations:
        {% for idx in annual_demand.index %}
            {{ idx }}.techs:
                biofuel_to_liquids:
                hydrogen_to_liquids:
                biofuel_to_diesel:
                biofuel_to_gas_and_liquids:
                biofuel_to_methanol:
                hydrogen_to_methanol:
                biofuel_to_methane:
                hydrogen_to_methane:
                electrolysis:
                dac:
                {% for carrier in carriers %}
                {% if annual_demand.loc[idx, carrier] > 1e-6 %}
                demand_industry_{{ carrier }}:
                {% else %}
                demand_industry_{{ carrier }}.exists: false
                {% endif %}
                {% endfor %}
        {% endfor %}
    new_hydrogen_storage:
        techs:
            hydrogen_storage: # based on [Danish energy agency, energy storage, 151a Hydrogen Storage - Tanks, 2050]
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
                    energy_cap: 0  # negate the costs given in the base euro-calliope definition (anaerobic digestion gas turbine)
                    om_prod: 0  # negate the costs given in the base euro-calliope definition (anaerobic digection gas turbine)
                    om_annual: 0  # negate the costs given in the base euro-calliope definition (anaerobic digection gas turbine)
                    om_con: {{ biofuel_fuel_cost * scaling_factors.specific_costs }} # {{ (1 / scaling_factors.specific_costs) | unit("EUR/MWh") }}
scenarios:
    industry_fuel_isolated: [annual_fuel_demand_isolated, industry_techs, new_hydrogen_storage, new_biofuel_supply, biofuel_maximum]
    industry_fuel_shared: [annual_fuel_demand_shared, industry_techs, new_hydrogen_storage, new_biofuel_supply, biofuel_maximum]
"""


def fill_constraint_template(
    path_to_annual_demand, path_to_result, model_year, scaling_factors,
    path_to_biofuel_costs, model_time, industry_carriers
):
    """Generate a file that represents links in Calliope."""
    annual_demand = util.read_tdf(path_to_annual_demand)
    annual_demand = (
        annual_demand
        .xs(('industry_demand', model_year), level=('dataset', 'year'))
        .sum(level=['id', 'end_use'])
        .unstack('end_use')
        .apply(util.filter_small_values, rel_tol=1e-4)  # remove excessively small values
        .loc[:, industry_carriers]
    )

    scaling_factors["specific_costs"] = scaling_factors["monetary"] / scaling_factors["power"]
    model_timedelta = util.get_timedelta(model_time, model_year)

    with open(path_to_biofuel_costs, "r") as f_biofuel_costs:
        biofuel_fuel_cost = float(f_biofuel_costs.readline())

    env = jinja2.Environment(lstrip_blocks=True, trim_blocks=True)
    env.filters["unit"] = filters.unit

    fuel = env.from_string(TEMPLATE).render(
        annual_demand=annual_demand,
        scaling_factors=scaling_factors,
        biofuel_fuel_cost=biofuel_fuel_cost,
        timedelta=model_timedelta,
        carriers=industry_carriers
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(fuel)


if __name__ == "__main__":
    fill_constraint_template(
        path_to_annual_demand=snakemake.input.annual_demand,
        path_to_result=snakemake.output[0],
        model_year=snakemake.params.model_year,
        scaling_factors=snakemake.params.scaling_factors,
        industry_carriers=snakemake.params.industry_carriers,
        model_time=snakemake.params.model_time,
        path_to_biofuel_costs=snakemake.input.biofuel_cost
    )
