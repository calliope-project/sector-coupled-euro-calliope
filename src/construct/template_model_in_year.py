from pathlib import Path

import jinja2


def parameterise_template(path_to_template, year, subset_time, resolution, path_to_result):
    """Copies template files, applying config parameters where relevant."""

    path_to_template = Path(path_to_template)
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(path_to_template.parent), lstrip_blocks=True, trim_blocks=True)
    assert all("year" in i for i in subset_time)
    subset_time = [i.replace("year", str(year)) for i in subset_time]
    rendered = env.get_template(path_to_template.name).render(
        subset_time=subset_time,
        year=year,
        resolution=resolution
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(rendered)


if __name__ == "__main__":
    parameterise_template(
        path_to_template=snakemake.input.template,
        year=snakemake.wildcards.year,
        subset_time=snakemake.params["subset_time"],
        resolution=snakemake.wildcards.resolution,
        path_to_result=snakemake.output[0]
    )
