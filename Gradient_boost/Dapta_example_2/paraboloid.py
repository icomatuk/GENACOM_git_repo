import os
from datetime import datetime
from pathlib import Path
import json
import shutil
from copy import deepcopy

HOSTNAME = os.getenv("HOSTNAME")

def compute(msg):
    
    resp = deepcopy(SETUP_DATA)
    inputs = resp.pop("input_data")
    
    c_inputs = {
            k: v
            for k, v in msg["inputs"]["design"].items()
        }
    
    x = c_inputs["x"][0]
    y = c_inputs["y"][0]

    outputs = resp.pop("output_data")
    outputs["design"]["f_xy"] = (x - 3.0) ** 2 + x * y + (y + 4.0) ** 2 - 3.0

    response = {}
    response["outputs"] = outputs

    options = resp.pop("options")  
    if options["get_grads"]:

        partials = resp.pop("partials")
        partials["f_xy"]["x"]["val"] = [2 * (x - 3.0) + y]
        partials["f_xy"]["y"]["val"] = [x + 2 * (y + 4.0)]
        response["partials"] = partials

    message = f"{datetime.now().strftime('%Y%m%d-%H%M%S')}: Compute paraboloid f(x:{str(x)},y:{str(y)}) = {str(outputs['design']['f_xy'])} on host {HOSTNAME}"
    resp["message"] = message

    return response

def setup(msg):
    """Editable setup function."""
    # import component parameters

    from design_optimisation import JSON_INPUT_PATH

    with open(JSON_INPUT_PATH, "r") as f:
        dapta_inputs = json.load(f)
        PARABOLOID_IN = [
            d for d in dapta_inputs["components"] if d["name"] == "paraboloid"
        ][0]
        open_mdao_in = [
            d for d in dapta_inputs["components"] if d["name"] == "open-mdao"
        ][0]
        outputs_folder = open_mdao_in["parameters"]["outputs_folder"]

    resp = {
        "input_data": {"design": {}, "implicit": {}, "setup": {}},
        "output_data": {"design": {}, "implicit": {}, "setup": {}},
        "options": {},
        "parameters": {
            "user_input_files": [],
            "inputs_folder_path": Path(__file__).parent,
            "outputs_folder_path": Path(__file__).parent / outputs_folder,
        },
    }

    run_folder = resp["parameters"]["outputs_folder_path"]

    for k, v in PARABOLOID_IN["parameters"].items():
        resp["parameters"][k] = v

    for k, v in PARABOLOID_IN["inputs"].items():
        if not "material" in k:  # only optimise on float type variables
            resp["input_data"]["design"][k] = v

    for k, v in PARABOLOID_IN["outputs"].items():
        resp["output_data"]["design"][k] = v

    for k, v in PARABOLOID_IN["options"].items():
        resp["options"][k] = v

    # initialise partials - required for OpenMDAO gradient-based optimisation
    resp["partials"] = {
        "f_xy": {
            "x": {"val": [0.0], "method": "exact"},
            "y": {"val": [0.0], "method": "exact"},
        }
    }

    resp["message"] = f"{datetime.now().strftime('%Y%m%d-%H%M%S')}: Setup completed."

    global SETUP_DATA
    SETUP_DATA = resp
    #print(SETUP_DATA)

    return resp