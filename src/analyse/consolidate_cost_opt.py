from pathlib import Path

import pandas as pd

from friendly_calliope.io import write_dpkg


def consolidate_spores(
    paths_to_friendlies, initial_keywords, name, description, model_resolution,
    path_to_output_friendly
):

    cost_opt_metrics = [get_friendly_metrics(_path) for _path in paths_to_friendlies]

    metrics_grouped = {}
    for metric in cost_opt_metrics[0].keys():
        if "weather_year" in cost_opt_metrics[0][metric].index.names:
            metrics_grouped[metric] = pd.concat([
                cost_opt_metric[metric] for cost_opt_metric in cost_opt_metrics
            ])
        else:
            metrics_grouped[metric] = cost_opt_metrics[0][metric]

    _metrics_to_dpkg(
        metrics_grouped, path_to_output_friendly,
        initial_keywords, name, description, model_resolution
    )



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
    metrics, path_to_output, initial_keywords, name, description, model_resolution
):
    keywords = initial_keywords + [
        "calliope", "Euro-Calliope",
        f"resolution={model_resolution}H",
        "model_year=2050",
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
        model_resolution=snakemake.wildcards.model_resolution,
        path_to_output_friendly=snakemake.output.friendly_data,
    )
