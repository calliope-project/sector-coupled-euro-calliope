import sys
import glob
import os

import calliope


def rerun_model(dir_path):
    calliope.set_log_verbosity()
    path_to_most_recent_spore_results = sorted(glob.glob(dir_path + "spore_*.nc"))[-1]

    most_recent_spore_num = int(
        os.path.basename(path_to_most_recent_spore_results)
        .replace("spore_", "").replace(".nc", "")
    )
    print(f"Continuing with results from SPORE {most_recent_spore_num}")
    cost_op_model = calliope.read_netcdf(dir_path + "spore_0.nc")
    cost_op_model.run_config["spores_options"]["skip_cost_op"] = True
    if most_recent_spore_num == 0:
        cost_op_model._model_data["group_cost_max"].loc[{"costs": "monetary"}] = (
            cost_op_model._model_data.cost.sum().item() *
            (1 + cost_op_model.run_config["spores_options"]["slack"])
        )
        if "spores" not in cost_op_model._model_data.coords:
            cost_op_model._model_data = cost_op_model._model_data.assign_coords(spores=("spores", [most_recent_spore_num]))
    else:
        spore_result = calliope.read_netcdf(path_to_most_recent_spore_results)
        for var in ["group_cost_max", "cost_energy_cap"]:
            cost_op_model._model_data[var] = spore_result._model_data[var]
        if "spores" not in spore_result._model_data.coords:
            cost_op_model._model_data = cost_op_model._model_data.assign_coords(spores=("spores", [most_recent_spore_num]))
        else:
            cost_op_model._model_data.coords["spores"] = [most_recent_spore_num]


    #cost_op_model.run_config["spores_options"]["save_per_spore_path"] = dir_path + "spore_{}.nc"
    if "objective_cost_class" in cost_op_model._model_data.data_vars:
        cost_op_model._model_data = cost_op_model._model_data.drop_vars(["objective_cost_class"])

    cost_op_model._model_data.attrs.pop("objective_function_value", None)
    cost_op_model._model_data.attrs.pop("termination_condition", None)
    print(
        "skip_cost_op: ",
        calliope.AttrDict.from_yaml_string(cost_op_model._model_data.attrs["run_config"]).spores_options.skip_cost_op
    )
    cost_op_model.run(force_rerun=True)


if __name__ == "__main__":
    rerun_model(sys.argv[1])
