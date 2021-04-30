import sys

import jinja2
import numpy as np
import pandas as pd

import util

sys.path.append('euro-calliope/src/')
import filters

TEMPLATE = """

overrides:
    transport_demand_exists:
        techs:
            {% for tech in annual_vehicles.columns %}
            demand_{{ tech }}_transport.exists: true
            {% endfor %}
    transport_ev_exists:
        techs:
            {% for tech in annual_vehicles.columns %}
            {{ tech }}_transport_ev.exists: true
            {% endfor %}
    transport_ice_exists:
        techs:
            {% for tech in annual_vehicles.columns %}
            {{ tech }}_transport_ice.exists: true
            {% endfor %}

    all_transport:
        locations:
        {% for idx in annual_vehicles.index %}
            {{ idx }}.techs:
                {% for tech in annual_vehicles.columns %}
                {{ tech }}_transport_ice:
                    constraints:
                        energy_eff: {{ efficiency_ice.loc[idx, tech] * scaling_factors.transport_efficiency }}  # {{ (1 / scaling_factors.transport_efficiency) | unit("Miokm/MWh") }}
                {% if ev_cap[tech].sum() > 0 %}
                {{ tech }}_transport_ev:
                    constraints:
                        energy_eff: {{ efficiency_ev.loc[idx, tech] * scaling_factors.transport_efficiency }}  # {{ (1 / scaling_factors.transport_efficiency) | unit("Miokm/MWh") }}
                        energy_cap_max: {{ ev_cap.loc[idx, tech] * scaling_factors.transport }}  # {{ (1 / scaling_factors.transport) | unit("Miokm") }}
                {% endif %}
                demand_{{ tech }}_transport:
                {% endfor %}
            {% endfor %}

    annual_transport_distance:
        group_constraints:
            {% for idx in annual_distance.index %}
            {% for tech in annual_distance.columns %}
            {{ tech }}_{{ idx }}_annual_distance:
                techs: [demand_{{ tech }}_transport]
                locs: [{{ idx }}]
                carrier_con_equals:
                    {{ tech }}_transport: {{ -1 * annual_distance.loc[idx, tech] * timedelta }}  # {{ (1 / scaling_factors.transport) | unit("Miokm") }}
            {% endfor %}

            {% endfor %}

    monthly_transport_demand_range:
        techs:
            {% for tech in annual_distance.columns %}
            {{ tech }}_transport_ev:
                constraints:
                    carrier_prod_per_month_min.{{ tech }}_transport: file=demand-min-{{ tech }}-ev.csv
                    carrier_prod_per_month_max.{{ tech }}_transport: file=demand-max-{{ tech }}-ev.csv
            {% endfor %}

    monthly_transport_demand_equality:
        techs:
            {% for tech in annual_distance.columns %}
            {{ tech }}_transport_ev:
                constraints:
                    carrier_prod_per_month_equals.{{ tech }}_transport: file=demand-equals-{{ tech }}-ev.csv
            {% endfor %}

scenarios:
    transport: [all_transport, annual_transport_distance, transport_demand_exists, transport_ev_exists, transport_ice_exists, monthly_transport_demand_equality]
"""


def fill_constraint_template(path_to_annual_demand, path_to_result, model_year, transport_params, scaling_factors, model_time):
    """Generate a file that represents links in Calliope."""
    annual_demand = util.read_tdf(path_to_annual_demand)

    def _get_data(annual_demand, dataset):
        _df = (
            annual_demand
            .drop(annual_demand.filter(regex='electricity').index)
            .xs((dataset, 'road', model_year), level=('dataset', 'cat_name', 'year'))
        )

        return _df.droplevel('unit').unstack('end_use')
    annual_distance = _get_data(annual_demand, 'transport_demand')
    annual_vehicles = _get_data(annual_demand, 'transport_vehicles')

    def _get_efficiency(x, fuel):
        return pd.Series(
            data=[1 / transport_params['efficiency'][i].get(fuel, np.nan) for i in x.index],
            index=x.index
        )

    def _wavg(x, weights):
        return pd.Series(
            data=np.average(x, weights=weights[x.columns], axis=1),
            index=x.index
        )

    efficiency_ice = annual_distance.apply(_get_efficiency, fuel='diesel', axis=1)
    efficiency_ev = annual_distance.apply(_get_efficiency, fuel='electricity', axis=1)
    ev_cap = annual_vehicles.mul(transport_params['ev_battery_capacity']).div({k: v.get('electricity', np.nan) for k, v in transport_params['efficiency'].items()}).fillna(0)

    if transport_params.get('group', None) is not None:
        # I.e. we group vehicle classes into either 'light' or 'heavy'
        efficiency_ice = efficiency_ice.groupby(transport_params['group'], axis=1).apply(_wavg, weights=annual_distance)
        efficiency_ev = efficiency_ev.groupby(transport_params['group'], axis=1).apply(_wavg, weights=annual_distance)
        ev_cap = ev_cap.groupby(transport_params['group'], axis=1).sum()
        annual_distance = annual_distance.groupby(transport_params['group'], axis=1).sum()
        annual_vehicles = annual_vehicles.groupby(transport_params['group'], axis=1).sum()

    model_timedelta = util.get_timedelta(model_time, model_year)
    scaling_factors["transport_efficiency"] = scaling_factors["transport"] / scaling_factors["power"]

    env = jinja2.Environment(lstrip_blocks=True, trim_blocks=True)
    env.filters["unit"] = filters.unit

    road_transport_constraints = env.from_string(TEMPLATE).render(
        annual_distance=annual_distance,
        annual_vehicles=annual_vehicles,
        efficiency_ice=efficiency_ice,
        efficiency_ev=efficiency_ev,
        ev_cap=ev_cap,
        scaling_factors=scaling_factors,
        timedelta=model_timedelta
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(road_transport_constraints)


if __name__ == "__main__":
    fill_constraint_template(
        path_to_annual_demand=snakemake.input.annual_demand,
        path_to_result=snakemake.output[0],
        model_year=snakemake.params.model_year,
        transport_params=snakemake.params.transport,
        model_time=snakemake.params.model_time,
        scaling_factors=snakemake.params.scaling_factors,
    )
