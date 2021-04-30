import argparse

import numpy as np

import calliope

EXCL_SCORE_COST = 100


def rerun_model(dir_path, n_spores, techs=None, tech_groups=None):
    calliope.set_log_verbosity()
    cost_op_model = calliope.read_netcdf(dir_path + "spore_0.nc")
    print("Loaded cost optimal model")
    cost_op_model.run_config["spores_options"]["skip_cost_op"] = True
    cost_op_model.run_config["spores_options"]["objective_cost_class"] = {"spores_score": 1, "excl_score": 10, "monetary": 0}
    cost_op_model.run_config["spores_options"]["spores_number"] = n_spores
    cost_op_model.run_config["spores_options"]["save_per_spore_path"] = dir_path + "spore_excl_{}_{{}}.nc".format([i for i in [techs, tech_groups] if i is not None])

    cost_op_model._model_data["group_cost_max"].loc[{"costs": "monetary"}] = (
        cost_op_model._model_data.cost.sum().item() *
        (1 + cost_op_model.run_config["spores_options"]["slack"])
    )
    if "spores" not in cost_op_model._model_data.coords:
        cost_op_model._model_data = cost_op_model._model_data.assign_coords(spores=("spores", [0]))

    if tech_groups is not None and techs is None:
        tech_groups = tech_groups.split(",")
        techs = np.unique(np.concatenate([
            cost_op_model._model_data.techs.loc[cost_op_model._model_data.inheritance.str.find(tech_group) > -1].values
            for tech_group in tech_groups
        ]))
    elif techs is not None:
        techs = techs.split(",")

    loc_techs = [
        f"{loc}::{tech}"
        for loc in cost_op_model._model_data.locs.values
        for tech in techs
    ]
    valid_loc_techs = cost_op_model._model_data.loc_techs_investment_cost.to_index().intersection(loc_techs)
    cost_op_model._model_data.cost_energy_cap.loc[{"costs": "excl_score", "loc_techs_investment_cost": valid_loc_techs}] = EXCL_SCORE_COST
    print(f"Running with exclusion scores activated for techs {techs} and a total of {len(loc_techs)} loc::techs")
    if "objective_cost_class" in cost_op_model._model_data.data_vars:
        cost_op_model._model_data = cost_op_model._model_data.drop_vars(["objective_cost_class"])

    cost_op_model._model_data.attrs.pop("objective_function_value", None)
    cost_op_model._model_data.attrs.pop("termination_condition", None)

    cost_op_model.run(force_rerun=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("spore_dir", help="Directory to load cost optimal model and save new spores")
    parser.add_argument("--n_spores", help="Number of exclusion spores to run", default=10, type=int)
    parser.add_argument("--techs", help="Tech to exclude")
    parser.add_argument("--tech_groups", help="Tech group to exclude")
    args = parser.parse_args()

    rerun_model(
        args.spore_dir, args.n_spores,
        techs=args.techs, tech_groups=args.tech_groups
    )
