import argparse
import glob
import os

import xarray as xr

import calliope


def rerun_model(dir_path, slack):
    calliope.set_log_verbosity()
    all_spores_results = glob.glob(dir_path + "spore_*")
    path_to_most_recent_spore_results = sorted([i for i in all_spores_results if "spore_excl-" not in i])[-1]

    most_recent_spore_num = int(
        os.path.basename(path_to_most_recent_spore_results)
        .replace("spore_", "").replace(".nc", "")
    )
    print(f"Continuing with results from SPORE {most_recent_spore_num}")

    cost_op_model = _prep_cost_op_model(dir_path, slack)
    relevant_loc_techs = _get_relevant_loc_techs(cost_op_model)
    if most_recent_spore_num == 0:
        new_scores = _get_new_scores(cost_op_model._model_data, cost_op_model)
    else:
        spore_result = calliope.read_netcdf(path_to_most_recent_spore_results)
        new_scores = _get_new_scores(spore_result._model_data, cost_op_model)
        cost_op_model._model_data["cost_energy_cap"].loc[{"costs": "spores_score"}] += spore_result._model_data["cost_energy_cap"].loc[{"costs": "spores_score"}]

    if "spores" not in cost_op_model._model_data.coords:
        cost_op_model._model_data = cost_op_model._model_data.assign_coords(spores=("spores", [most_recent_spore_num]))
    else:
        cost_op_model._model_data.coords["spores"] = [most_recent_spore_num]

    if "objective_cost_class" in cost_op_model._model_data.data_vars:
        cost_op_model._model_data = cost_op_model._model_data.drop_vars(["objective_cost_class"])

    print(f"SPORES scores starting out summing to {cost_op_model._model_data.cost_energy_cap.loc[{'costs': 'spores_score'}].sum().item()}")
    cost_op_model._model_data["cost_energy_cap"].loc[{"costs": "spores_score", "loc_techs_investment_cost": relevant_loc_techs}] += new_scores
    print(f"SPORES scores being sent to optimisation summing to {cost_op_model._model_data.cost_energy_cap.loc[{'costs': 'spores_score'}].sum().item()}")

    assert len(new_scores) == len(relevant_loc_techs)
    assert len((cost_op_model._model_data["cost_energy_cap"].loc[{"costs": "spores_score"}].dropna("loc_techs_investment_cost"))) == len(relevant_loc_techs)
    print(f"Running SPORES with slack of {cost_op_model.run_config['spores_options']['slack']}")
    cost_op_model.run(force_rerun=True)


def _get_new_scores(model_data, cost_op_model):
    relevant_loc_techs = _get_relevant_loc_techs(cost_op_model)

    # First, we ensure that we don't erroneously score techs that have to be there no matter what.
    if "energy_cap_min" in cost_op_model._model_data.data_vars.keys():
        min_cap = cost_op_model._model_data.energy_cap_min.fillna(0)
    else:
        min_cap = 0

    if "energy_cap_equals" in cost_op_model._model_data.data_vars.keys():
        forced_energy_cap = cost_op_model._model_data.energy_cap_equals.fillna(min_cap)
        print(forced_energy_cap.to_series().filter(regex="transmission"))
    else:
        forced_energy_cap = min_cap

    return xr.where(model_data.energy_cap - forced_energy_cap > 1e-3, 100, 0).loc[relevant_loc_techs].drop_vars("loc_techs")


def _get_relevant_loc_techs(cost_op_model):
    return (
        cost_op_model
        ._model_data
        .cost_energy_cap
        .loc[{"costs": "spores_score"}]
        .dropna("loc_techs_investment_cost")
        .loc_techs_investment_cost
    )


def _prep_cost_op_model(dir_path, slack):
    cost_op_model = calliope.read_netcdf(dir_path + "spore_0.nc")
    cost_op_model.run_config["spores_options"]["save_per_spore_path"] = dir_path + "spore_{}.nc"
    cost_op_model.run_config["spores_options"]["skip_cost_op"] = True
    cost_op_model.run_config["spores_options"]["slack"] = slack
    cost_op_model.run_config["objective_options"]["cost_class"] = {"spores_score": 1, "excl_score": 0, "monetary": 0}
    cost_op_model._model_data["group_cost_max"].loc[{"costs": "monetary"}] = (
        cost_op_model._model_data.cost.sum().item() *
        (1 + cost_op_model.run_config["spores_options"]["slack"])
    )
    cost_op_model._model_data.cost_depreciation_rate.loc[{"costs": "spores_score"}] = (
        cost_op_model._model_data.cost_depreciation_rate.loc[{"costs": "spores_score"}]
        .clip(min=1, max=1)
    )

    cost_op_model._model_data.attrs.pop("objective_function_value", None)
    cost_op_model._model_data.attrs.pop("termination_condition", None)
    return cost_op_model


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("spore_dir", help="Directory to load cost optimal model and save new spores")
    parser.add_argument("--slack", help="SPORES slack", default=0.1, type=float)
    args = parser.parse_args()
    rerun_model(args.spore_dir, args.slack)
