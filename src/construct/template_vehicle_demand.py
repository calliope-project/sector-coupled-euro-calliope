import jinja2

import util

TEMPLATE = """

overrides:
    all_transport:
        locations:
        {% for idx in annual_vehicles.index %}
            {{ idx }}.techs:
                hdv_transport_ice:
                ldv_transport_ice:
                ldv_transport_ev:
                bus_transport_ice:
                bus_transport_ev:
                motorcycle_transport_ice:
                motorcycle_transport_ev:
                passenger_car_transport_ice:
                passenger_car_transport_ev:
                demand_passenger_car_transport:
                demand_bus_transport:
                demand_motorcycle_transport:
                demand_ldv_transport:
                demand_hdv_transport:
            {% endfor %}

    max_ev_battery_charge_capacity:
        locations:
        {% for idx in annual_vehicles.index %}
            {{ idx }}.techs:
            {% for tech in annual_vehicles.columns %}
                {% if ev_battery_capacity[tech] is defined and efficiency[tech].electricity is defined %}
                {{ tech }}_transport_ev.constraints.energy_cap_max: {{ annual_vehicles.loc[idx, tech] * ev_battery_capacity[tech] / (efficiency[tech].electricity / scaling) }}  # {{ distance_unit }}
                {% endif %}
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
                    {{ tech }}_transport: {{ -1 * annual_distance.loc[idx, tech] * timedelta }}  # {{ distance_unit }}
            {% endfor %}

            {% endfor %}
scenarios:
    transport: [all_transport, max_ev_battery_charge_capacity, annual_transport_distance]
"""


def fill_constraint_template(path_to_annual_demand, path_to_result, model_year, transport_params, scaling, model_time):
    """Generate a file that represents links in Calliope."""
    annual_demand = util.read_tdf(path_to_annual_demand)

    def _get_data(annual_demand, dataset):
        _df = (
            annual_demand
            .drop(annual_demand.filter(regex='electricity').index)
            .xs((dataset, 'road', model_year), level=('dataset', 'cat_name', 'year'))
        )
        _unit = _df.index.get_level_values('unit').unique()
        assert len(_unit) == 1
        return _df.droplevel('unit').unstack('end_use'), _unit[0]
    annual_distance, distance_unit = _get_data(annual_demand, 'transport_demand')
    annual_vehicles, vehicles_unit = _get_data(annual_demand, 'transport_vehicles')
    model_timedelta = util.get_timedelta(model_time, model_year)

    template = jinja2.Environment(lstrip_blocks=True, trim_blocks=True).from_string(TEMPLATE)
    road_transport_constraints = template.render(
        annual_distance=annual_distance,
        distance_unit=distance_unit,
        annual_vehicles=annual_vehicles,
        ev_battery_capacity=transport_params['ev_battery_capacity'],
        efficiency=transport_params['efficiency'],
        scaling=scaling,  # to scale the distance part of efficiency
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
        scaling=snakemake.params.scaling,
    )
