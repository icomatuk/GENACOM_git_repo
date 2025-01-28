import shutil
from pathlib import Path
import json
from datetime import datetime
from copy import deepcopy

TEST_MSG_SETUP = {"component": "link_variables"}
TEST_MSG_COMPUTE = {
    "component": "link_variables",
    "inputs": {"design": {}},
    "get_grads": False,
    "get_outputs": True,
}
TEST_JSON_INPUT_PATH = None

SETUP_DATA = {}


def setup(msg):
    """Editable setup function."""
    # import component parameters
    if TEST_JSON_INPUT_PATH:
        JSON_INPUT_PATH = TEST_JSON_INPUT_PATH
    else:
        from design_optimisation import JSON_INPUT_PATH

    with open(JSON_INPUT_PATH, "r") as f:
        dapta_inputs = json.load(f)
        LINK_VARIABLES_IN = [
            d for d in dapta_inputs["components"] if d["name"] == "link_variables"
        ][0]
        open_mdao_in = [
            d for d in dapta_inputs["components"] if d["name"] == "open-mdao"
        ][0]
        outputs_folder = open_mdao_in["parameters"]["outputs_folder"]

    resp = {
        "input_data": {"design": {}, "implicit": {}, "setup": {}},
        "output_data": {"design": {}, "implicit": {}, "setup": {}},
        "parameters": {
            "user_input_files": [],
            "inputs_folder_path": Path(__file__).parent,
            "outputs_folder_path": Path(__file__).parent / outputs_folder,
        },
    }

    run_folder = resp["parameters"]["outputs_folder_path"]

    for k, v in LINK_VARIABLES_IN["parameters"].items():
        resp["parameters"][k] = v

    for k, v in LINK_VARIABLES_IN["inputs"].items():
        if not "material" in k:  # only optimise on float type variables
            resp["input_data"]["design"][k] = v

    for k, v in LINK_VARIABLES_IN["outputs"].items():
        resp["output_data"]["design"][k] = v

    resp["message"] = f"{datetime.now().strftime('%Y%m%d-%H%M%S')}: Setup completed."

    global SETUP_DATA
    SETUP_DATA = resp
    # print(SETUP_DATA)

    return resp


def compute(msg):
    resp = deepcopy(SETUP_DATA)
    run_folder = resp["parameters"]["outputs_folder_path"]
    resp.pop("input_data")

    # read inputs from object
    composite_inputs_parsed = {}
    c_inputs = {
        k: v
        for k, v in msg["inputs"]["design"].items()
        if k.startswith("REF_composite.")
    }
    for k, v in c_inputs.items():
        if isinstance(v, list):
            v = v[0]
        kchain = k.split(".")
        if len(kchain) == 5 and kchain[4] in ["thickness", "angle"]:
            pass
        elif len(kchain) == 6 and kchain[5] == "inc_angle":
            pass               
        else:
            raise ValueError(f"Input variable of type {'.'.join(kchain[2:])} is not recognised.")
        instance = kchain[1]
        if not instance in composite_inputs_parsed:
            composite_inputs_parsed[instance] = {}     
        composite_inputs_parsed[instance][".".join(kchain[2:])] = float(v) 

    # define outputs
    outputs = resp.pop("output_data")
    c_outputs = {
        k: v for k, v in outputs["design"].items() if k.startswith("composite.")
    }
    for output in c_outputs:
        kchain = output.split(".")
        if len(kchain) == 5 and kchain[4] in ["thickness", "angle"]:
            instance = kchain[1]
            layup = kchain[2]
            ply = int(kchain[3])
            v_type = kchain[4]
            formula = resp["parameters"]["composite"][instance][layup]["plies"][str(ply)][v_type]
        elif len(kchain) == 6 and kchain[5] == "inc_angle":
            instance = kchain[1]
            layup = kchain[2]
            ply = kchain[3]
            control_point = int(kchain[4])
            v_type = kchain[5]
            formula = resp["parameters"]["rts_plies"][instance][ply][str(control_point)][v_type]
        else:
            raise ValueError(f"Output variable of type {'.'.join(kchain[2:])} is not recognised.")
        
        # lookup link to variables
        outputs["design"][output] = get_value_from_formula(
            formula, composite_inputs_parsed[instance]
        )

    response = {}
    response["outputs"] = outputs
    resp["message"] = f"{datetime.now().strftime('%Y%m%d-%H%M%S')}: Compute completed."

    return response


def get_value_from_formula(formula, inputs):
    if formula.startswith("equal("):
        fp = formula.split("(")[1]
        fp = fp.replace(")", "")
        if "," in fp:
            raise ValueError(
                f"Equal formula should only have one input parameter: {fp}"
            )
        return inputs.get(fp)
    
    if formula.startswith("negative("):
        fp = formula.split("(")[1]
        fp = fp.replace(")", "")
        if "," in fp:
            raise ValueError(
                f"Equal formula should only have one input parameter: {fp}"
            )
        return -inputs.get(fp)

    if formula.startswith("linear("):
        fp = formula.split("(")[1]
        fp = fp.replace(")", "")
        fp = fp.split(",")
        if len(fp) != 3:
            raise ValueError(f"Linear formula should have three input parameters: {fp}")

        v1 = inputs.get(fp[0])
        v2 = inputs.get(fp[1])
        ratio = float(fp[2])
        if not ratio >= 0 or not ratio <= 1:
            raise ValueError(f"ratio value should be float between 0. and 1.")
        v = (v2 - v1) * ratio + v1

        return v      

    raise ValueError(f"Cannot understand formula {f}")


if __name__ == "__main__":
    json_input_file = "Study_3_4_c.json"

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
    TEST_JSON_INPUT_PATH = Path(outputs_folder) / "dapta_inputs.json"

    resp = setup(TEST_MSG_SETUP)
    print(resp)

    TEST_MSG_COMPUTE["inputs"]["design"] = SETUP_DATA["input_data"]["design"]
    resp = compute(TEST_MSG_COMPUTE)
    print(resp)

    print("TEST linked_variables.py passed.")
