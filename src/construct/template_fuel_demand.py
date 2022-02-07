import jinja2

import util

from eurocalliopelib import filters

TEMPLATE = """
overrides:
    annual_fuel_demand:
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
                syn_diesel_converter:
                syn_methane_converter:
                syn_kerosene_converter:
                syn_methanol_converter:
                {% for carrier in carriers %}
                {% if annual_demand.loc[idx, carrier] > 1e-6 %}
                demand_industry_{{ carrier }}:
                {% else %}
                demand_industry_{{ carrier }}.exists: false
                {% endif %}
                {% endfor %}
        {% endfor %}
scenarios:
    industry_fuel: [annual_fuel_demand, industry_techs]
"""


def fill_constraint_template(
    path_to_annual_demand, path_to_result, year, scaling_factors,
    model_time, industry_carriers
):
    """Generate a file that represents links in Calliope."""
    annual_demand = util.read_tdf(path_to_annual_demand)
    annual_demand = (
        annual_demand.xs(('industry_demand', int(year)), level=('dataset', 'year')).sum(level=['id', 'end_use']).unstack('end_use')
        .reindex(industry_carriers, axis=1)
        .fillna(0)
        .apply(util.filter_small_values, rel_tol=1e-4)  # remove excessively small values
    )

    scaling_factors["specific_costs"] = scaling_factors["monetary"] / scaling_factors["power"]
    model_timedelta = util.get_timedelta(model_time, year)

    env = jinja2.Environment(lstrip_blocks=True, trim_blocks=True)
    env.filters["unit"] = filters.unit

    fuel = env.from_string(TEMPLATE).render(
        annual_demand=annual_demand,
        scaling_factors=scaling_factors,
        timedelta=model_timedelta,
        carriers=industry_carriers
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(fuel)


if __name__ == "__main__":
    fill_constraint_template(
        path_to_annual_demand=snakemake.input.annual_demand,
        path_to_result=snakemake.output[0],
        year=snakemake.wildcards.year,
        scaling_factors=snakemake.params.scaling_factors,
        industry_carriers=snakemake.params.industry_carriers,
        model_time=snakemake.params.model_time
    )
