from pathlib import Path

import jinja2


def parameterise_template(path_to_template, shares, subset_time, path_to_result):
    """Copies template files, applying config parameters where relevant."""

    path_to_template = Path(path_to_template)
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(path_to_template.parent), lstrip_blocks=True, trim_blocks=True)
    env.globals["int"] = int

    rendered = env.get_template(path_to_template.name).render(
        shares=shares,
        subset_time=subset_time
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(rendered)


if __name__ == "__main__":
    parameterise_template(
        path_to_template=snakemake.input.template,
        shares=snakemake.params["shares"],
        subset_time=snakemake.params["subset_time"],
        path_to_result=snakemake.output[0]
    )
