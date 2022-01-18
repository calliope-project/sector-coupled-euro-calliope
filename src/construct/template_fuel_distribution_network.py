import pandas as pd

import jinja2

TEMPLATE = """
overrides:
    synfuel_transmission:
        techs:
            syn_diesel_distribution_export:
                essentials:
                    carrier: syn_diesel
                    parent: demand
                constraints:
                    resource: -.inf
                    force_resource: false

            syn_methane_distribution_export:
                essentials:
                    carrier: syn_methane
                    parent: demand
                constraints:
                    resource: -.inf
                    force_resource: false

            syn_kerosene_distribution_export:
                essentials:
                    carrier: syn_kerosene
                    parent: demand
                constraints:
                    resource: -.inf
                    force_resource: false

            syn_methanol_distribution_export:
                essentials:
                    carrier: syn_methanol
                    parent: demand
                constraints:
                    resource: -.inf
                    force_resource: false

            syn_diesel_distribution_import:
                essentials:
                    carrier: syn_diesel
                    parent: supply
                constraints:
                    resource: .inf

            syn_methane_distribution_import:
                essentials:
                    carrier: syn_methane
                    parent: supply
                constraints:
                    resource: .inf

            syn_kerosene_distribution_import:
                essentials:
                    carrier: syn_kerosene
                    parent: supply
                constraints:
                    resource: .inf

            syn_methanol_distribution_import:
                essentials:
                    carrier: syn_methanol
                    parent: supply
                constraints:
                    resource: .inf

        locations:
            {{ regions|join(",") }}:
                techs:
                    syn_diesel_distribution_export:
                    syn_methane_distribution_export:
                    syn_diesel_distribution_import:
                    syn_methane_distribution_import:
                    syn_kerosene_distribution_export:
                    syn_methanol_distribution_export:
                    syn_kerosene_distribution_import:
                    syn_methanol_distribution_import:


"""


def parametrise_template(path_to_regions, path_to_output):
    regions = pd.read_csv(path_to_regions, index_col=0).index
    env = jinja2.Environment(lstrip_blocks=True, trim_blocks=True)

    rendered = env.from_string(TEMPLATE).render(
        regions=regions
    )
    with open(path_to_output, "w") as result_file:
        result_file.write(rendered)


if __name__ == "__main__":
    parametrise_template(
        path_to_regions=snakemake.input.regions,
        path_to_output=snakemake.output[0]
    )
