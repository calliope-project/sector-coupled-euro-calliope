import jinja2
import pandas as pd

from eurocalliopelib import filters

TEMPLATE = """
links:
    {% for row_id, row in gtcs.iterrows() %}
    {{ row_id[0] }},{{ row_id[1] }}.techs:
        {% for link, data in row.iteritems() %}
        {% if data > 0 %}
        {{ link[0] }}_{{ link[1].lower().replace('-', '_') }}_transmission:
            constraints:
                energy_cap_min: {{ data * scaling_factors.power }}  # {{ (1 / scaling_factors.power) | unit("MW") }}
            costs.monetary.energy_cap: {{ costs["{}-{}".format(link[0], link[1].lower())] * scaling_factors.specific_costs }}  # {{ (1 / scaling_factors.specific_costs) | unit("EUR/MW") }}
        {% if link[2] == "one-way" %}
            one_way: true
        {% endif %}
        {% endif %}
        {% endfor %}
    {% endfor %}

overrides:
    transmission_om_cost:
        tech_groups:
            dc_transmission.costs.monetary.om_prod: 1e-5  # dummy cost
            ac_transmission.costs.monetary.om_prod: 1e-5  # dummy cost

    link_cap_20:
        tech_groups:
            dc_transmission.constraints.energy_cap_max: 0.2
            ac_transmission.constraints.energy_cap_max: 0.2

    link_cap_40:
        tech_groups:
            dc_transmission.constraints.energy_cap_max: 0.4
            ac_transmission.constraints.energy_cap_max: 0.4

    link_cap_10_40:
        tech_groups:
            dc_transmission.constraints.energy_cap_max: 0.1
            ac_transmission.constraints.energy_cap_max: 0.4

    link_cap_1x:
        links:
            {% for row_id, row in gtcs.iterrows() %}
            {{ row_id[0] }},{{ row_id[1] }}.techs:
                {% for link, data in row.iteritems() %}
                {% if data > 0 %}
                {{ link[0] }}_{{ link[1].lower().replace('-', '_') }}_transmission.constraints.energy_cap_equals: {{ data * scaling_factors.power }}  # {{ (1 / scaling_factors.power) | unit("MW") }}
                {% endif %}
                {% endfor %}
            {% endfor %}

    link_cap_2x:
        links:
            {% for row_id, row in gtcs.iterrows() %}
            {{ row_id[0] }},{{ row_id[1] }}.techs:
                {% for link, data in row.iteritems() %}
                {% if data > 0 %}
                {{ link[0] }}_{{ link[1].lower().replace('-', '_') }}_transmission.constraints.energy_cap_max: {{ data * 2 * scaling_factors.power }}  # {{ (1 / scaling_factors.power) | unit("MW") }}
                {% endif %}
                {% endfor %}
            {% endfor %}

    link_cap_5x:
        links:
            {% for row_id, row in gtcs.iterrows() %}
            {{ row_id[0] }},{{ row_id[1] }}.techs:
                {% for link, data in row.iteritems() %}
                {% if data > 0 %}
                {{ link[0] }}_{{ link[1].lower().replace('-', '_') }}_transmission.constraints.energy_cap_max: {{ data * 5 * scaling_factors.power }}  # {{ (1 / scaling_factors.power) | unit("MW") }}
                {% endif %}
                {% endfor %}
            {% endfor %}

    link_cap_10x:
        links:
            {% for row_id, row in gtcs.iterrows() %}
            {{ row_id[0] }},{{ row_id[1] }}.techs:
                {% for link, data in row.iteritems() %}
                {% if data > 0 %}
                {{ link[0] }}_{{ link[1].lower().replace('-', '_') }}_transmission.constraints.energy_cap_max: {{ data * 10 * scaling_factors.power }}  # {{ (1 / scaling_factors.power) | unit("MW") }}
                {% endif %}
                {% endfor %}
            {% endfor %}

    link_cap_50x:
        links:
            {% for row_id, row in gtcs.iterrows() %}
            {{ row_id[0] }},{{ row_id[1] }}.techs:
                {% for link, data in row.iteritems() %}
                {% if data > 0 %}
                {{ link[0] }}_{{ link[1].lower().replace('-', '_') }}_transmission.constraints.energy_cap_max: {{ data * 50 * scaling_factors.power }}  # {{ (1 / scaling_factors.power) | unit("MW") }}
                {% endif %}
                {% endfor %}
            {% endfor %}

    link_cap_100x:
        links:
            {% for row_id, row in gtcs.iterrows() %}
            {{ row_id[0] }},{{ row_id[1] }}.techs:
                {% for link, data in row.iteritems() %}
                {% if data > 0 %}
                {{ link[0] }}_{{ link[1].lower().replace('-', '_') }}_transmission.constraints.energy_cap_max: {{ data * 100 * scaling_factors.power }}  # {{ (1 / scaling_factors.power) | unit("MW") }}
                {% endif %}
                {% endfor %}
            {% endfor %}

    link_cap_500x:
        links:
            {% for row_id, row in gtcs.iterrows() %}
            {{ row_id[0] }},{{ row_id[1] }}.techs:
                {% for link, data in row.iteritems() %}
                {% if data > 0 %}
                {{ link[0] }}_{{ link[1].lower().replace('-', '_') }}_transmission.constraints.energy_cap_max: {{ data * 500 * scaling_factors.power }}  # {{ (1 / scaling_factors.power) | unit("MW") }}
                {% endif %}
                {% endfor %}
            {% endfor %}


    link_cap_dynamic:
        links:
            {% for row_id, row in gtcs.iterrows() %}
            {{ row_id[0] }},{{ row_id[1] }}.techs:
                {% for link, data in row.iteritems() %}
                {% if data > 0 and data <= 1000 %}
                {{ link[0] }}_{{ link[1].lower().replace('-', '_') }}_transmission.constraints.energy_cap_max: {{ data * 40 * scaling_factors.power }}  # {{ (1 / scaling_factors.power) | unit("MW") }}
                {% elif data > 1000 and data <= 5000 %}
                {{ link[0] }}_{{ link[1].lower().replace('-', '_') }}_transmission.constraints.energy_cap_max: {{ data * 10 * scaling_factors.power }}  # {{ (1 / scaling_factors.power) | unit("MW") }}
                {% elif data > 5000 and data <= 10000 %}
                {{ link[0] }}_{{ link[1].lower().replace('-', '_') }}_transmission.constraints.energy_cap_max: {{ data * 5 * scaling_factors.power }}  # {{ (1 / scaling_factors.power) | unit("MW") }}
                {% elif data > 10000 and data <= 15000 %}
                {{ link[0] }}_{{ link[1].lower().replace('-', '_') }}_transmission.constraints.energy_cap_max: {{ data * 3 * scaling_factors.power }}  # {{ (1 / scaling_factors.power) | unit("MW") }}
                {% elif data > 15000 %}
                {{ link[0] }}_{{ link[1].lower().replace('-', '_') }}_transmission.constraints.energy_cap_max: {{ data * 2 * scaling_factors.power }}  # {{ (1 / scaling_factors.power) | unit("MW") }}
                {% endif %}
                {% endfor %}
            {% endfor %}
"""


def generate_links(path_to_gtc, scaling_factors, path_to_result, resolution, costs):
    """Generate a file that represents links in Calliope."""
    gtcs = _read_gtcs(path_to_gtc, resolution)

    scaling_factors["specific_costs"] = scaling_factors["monetary"] / scaling_factors["power"]

    env = jinja2.Environment(lstrip_blocks=True, trim_blocks=True)
    env.filters["unit"] = filters.unit

    links = env.from_string(TEMPLATE).render(
        gtcs=gtcs,
        scaling_factors=scaling_factors,
        costs=costs
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(links)


def _read_gtcs(path_to_gtc, resolution):

    gtcs = (
        pd.read_csv(path_to_gtc, header=0)
        .set_index(['from', 'to', 'current', 'type', 'direction'])
        .drop('comment', axis=1)
        .squeeze()
        .unstack(['current', 'type', 'direction'])
    )
    if resolution == 'national':
        gtcs.index = gtcs.index.map(lambda x: tuple(i.split('_')[0] for i in x))
        # Remove internal links and sum links across the same borders
        gtcs = gtcs.groupby(level=[0, 1]).sum(min_count=1).drop([i for i in gtcs.index.values if i[0] == i[1]])
        gtc_duplicates = gtcs.reorder_levels([1, 0]).rename_axis(index=['from', 'to']).reindex(gtcs.index)
        gtcs = gtcs.fillna(0).add(gtc_duplicates.fillna()).rename(columns={"one-way": "two-way"}).groupby([0, 1, 2], axis=1).sum()
        to_drop = []
        for row_id in gtc_duplicates.index:
            if (row_id[1], row_id[0]) in gtc_duplicates.loc[:row_id].index:
                to_drop.append(row_id)
        gtcs = gtcs.drop(to_drop)

    return gtcs


if __name__ == "__main__":
    generate_links(
        path_to_gtc=snakemake.input.gtc,
        path_to_result=snakemake.output[0],
        scaling_factors=snakemake.params.scaling_factors,
        costs=snakemake.params.costs,
        resolution=snakemake.wildcards[0]
    )
