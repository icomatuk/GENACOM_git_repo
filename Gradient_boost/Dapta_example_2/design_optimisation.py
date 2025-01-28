""" Main script controlling the design optimisation task."""
import shutil
from pathlib import Path
from datetime import datetime
import json

JSON_INPUT_PATH = "script_unconstrained_fd.json"

from open_mdao import compute, post_process_optimisation 


def main(json_input_file):
    print(
        f"{datetime.now().strftime('%Y%m%d-%H%M%S')}: Starting design optimisation from {json_input_file} \n"
    )

    # import component parameters
    with open(json_input_file, "r") as f:
        dapta_inputs = json.load(f)
        open_mdao_in = [
            d for d in dapta_inputs["components"] if d["name"] == "open-mdao"
        ][0]
        outputs_folder = open_mdao_in["parameters"]["outputs_folder"]

    # copy input file into output folder
    shutil.copy2(
        json_input_file,
        Path(outputs_folder) / "dapta_inputs.json",
    )

    # define global variables
    resp = {
        "inputs": {"design": {}, "implicit": {}, "setup": {}},
        "outputs": {"design": {}, "implicit": {}, "setup": {}},
        "parameters": {
            "user_input_files": [],
            "inputs_folder_path": Path(__file__).parent,
            "outputs_folder_path": Path(__file__).parent / outputs_folder,
        },
    }

    for k, v in open_mdao_in["parameters"].items():
        resp["parameters"][k] = v

    # launch optimisation driver
    outputs = compute(**resp)


def plot_histories_from_file(json_input_file, recorder_file):
    """Re-plot results - same as switching on plot_history in the optimisation parameters."""

    # import component parameters
    with open(json_input_file, "r") as f:
        dapta_inputs = json.load(f)
        open_mdao_in = [
            d for d in dapta_inputs["components"] if d["name"] == "open-mdao"
        ][0]
        outputs_folder = open_mdao_in["parameters"]["outputs_folder"]

    # define global variables
    resp = {
        "inputs": {"design": {}, "implicit": {}, "setup": {}},
        "outputs": {"design": {}, "implicit": {}, "setup": {}},
        "parameters": {
            "user_input_files": [],
            "inputs_folder_path": Path(__file__).parent,
            "outputs_folder_path": Path(__file__).parent / outputs_folder,
        },
    }

    for k, v in open_mdao_in["parameters"].items():
        resp["parameters"][k] = v

    run_folder = resp["parameters"]["outputs_folder_path"]

    r_name = run_folder / recorder_file
    if not r_name.is_file():
        raise ValueError("can't find recorder file.")

    post_process_optimisation(resp["parameters"], run_folder, r_name)


if __name__ == "__main__":
    main(JSON_INPUT_PATH)
    # plot_histories_from_file(JSON_INPUT_PATH, "om_problem_recorder_20241015-101608.sqlite")
