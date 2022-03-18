from pathlib import Path

import pandas as pd

from friendly_calliope.io import write_dpkg


def consolidate_spores(
    paths_to_friendlies, initial_keywords, name, description, slack,
    path_to_output_friendly, path_to_processed_spores
):

    spores_metrics = {
        Path(spore): get_friendly_metrics(spore) for spore in paths_to_friendlies
    }
    processed_spores = pd.Series(list(spores_metrics.keys()))

    spores_metrics_grouped = {}
    first_spore = spores_metrics[list(spores_metrics.keys())[0]]
    for metric in first_spore.keys():
        if "spore" in first_spore[metric].index.names:
            spores_metrics_grouped[metric] = pd.concat([
                spore[metric].droplevel("spore")
                for spore in spores_metrics.values()
            ], keys=processed_spores.index, names=["spore"])
        else:
            spores_metrics_grouped[metric] = first_spore[metric]

    _metrics_to_dpkg(
        spores_metrics_grouped, path_to_output_friendly,
        initial_keywords, name, description, slack
    )

    processed_spores.to_csv(path_to_processed_spores)


def get_friendly_metrics(paths_to_friendly_data):
    metrics = {}
    metric_paths = Path(paths_to_friendly_data).glob("data/*.csv")
    for metric_path in metric_paths:
        metric_df = pd.read_csv(metric_path, index_col=False)
        metric_series = metric_df.set_index([i for i in metric_df.columns[:-1]]).squeeze(axis=1)
        assert isinstance(metric_series, pd.Series), f"{metric_path}"
        metrics[metric_series.name] = metric_series
    return metrics


def _metrics_to_dpkg(
    metrics, path_to_output, initial_keywords, name, description, slack
):
    if slack == "_5slack":
        cost_relaxation = 5
    elif slack == "_15slack":
        cost_relaxation = 15
    else:
        cost_relaxation = 10

    if slack == "_demand_update":
        initial_keywords.append("Projected 2050 demands")

    keywords = initial_keywords + [
        "calliope", "Euro-Calliope",
        "resolution=2H",
        "weather_year=2018",
        "model_year=2050",
        f"cost_relaxation={cost_relaxation}%"
    ]
    meta = {
        "name": name,
        "description": description,
        "keywords": keywords,
        "licenses": "CC-BY-4.0"
    }

    write_dpkg(metrics, path_to_output, meta, include_timeseries_data=True)


if __name__ == '__main__':
    consolidate_spores(
        paths_to_friendlies=snakemake.input.all_friendly_files,
        initial_keywords=snakemake.params.initial_keywords,
        name=snakemake.params.name,
        description=snakemake.params.description,
        slack=snakemake.wildcards.slack,
        path_to_output_friendly=snakemake.output.friendly_data,
        path_to_processed_spores=snakemake.output.processed_spores
    )
