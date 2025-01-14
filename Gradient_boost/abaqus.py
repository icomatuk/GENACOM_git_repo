import os
import subprocess
import shutil
from pathlib import Path
import json
from datetime import datetime
from copy import deepcopy
from math import sqrt, atan, degrees
from time import sleep
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import os
# from icomat_rts_cli import rts_analysis as rts
import rts_composite as rts

PLOT_OPTION = rts.PLOT_OPTION
LOCAL_EXECUTES = {
    "ABAQUS": "abaqus",
    "ABAQUS_POST": "abaqus python",
}
TEST_MSG_SETUP = {"component": "abaqus"}
TEST_MSG_COMPUTE = {
    "component": "abaqus", 
    "inputs": {"design": {}},
    "get_grads": False,
    "get_outputs": True,
}
TEST_JSON_INPUT_PATH = None

PLIES_EXCLUDED_FROM_THICKNESS_RATIO = rts.PLIES_EXCLUDED_FROM_THICKNESS_RATIO

SETUP_DATA = rts.SETUP_DATA

OUTPUT_type = {"MASS": np.float64,
               "U": np.float64,
               "T-RATIO-MIN": np.float64,
               "T-RATIO-MAX": np.float64,
               "TSAIW": np.float64,
               "E": np.float64
               }

def setup(msg):
    """Editable setup function."""

    # import component parameters
    if TEST_JSON_INPUT_PATH:
        JSON_INPUT_PATH = TEST_JSON_INPUT_PATH
    else:
        from design_optimisation import JSON_INPUT_PATH        # Theo revisit

    with open(JSON_INPUT_PATH, "r") as f:
        dapta_inputs = json.load(f)
        ABAQUS_IN = [d for d in dapta_inputs["components"] if d["name"] == "abaqus"][0]
        open_mdao_in = [
            d for d in dapta_inputs["components"] if d["name"] == "open-mdao"
        ][0]
        outputs_folder = open_mdao_in["parameters"]["outputs_folder"]

        # Theo
        # SETTINGS = dapta_inputs["settings"]    
        # rts.CPUs = SETTINGS["cpus"] 
        # print("CPUs used: {0}".format(rts.CPUs))

    resp = {
        "input_data": {"design": {}, "implicit": {}, "setup": {}},
        "output_data": {"design": {}, "implicit": {}, "setup": {}},
        "parameters": {
            "user_input_files": [],
            "inputs_folder_path": Path(os.getcwd()), #Path(__file__).parent,
            "outputs_folder_path": Path(os.getcwd()) / outputs_folder, #Path(__file__).parent / outputs_folder,
        },
    }

    run_folder = resp["parameters"]["outputs_folder_path"]

    for k, v in ABAQUS_IN["parameters"].items():
        resp["parameters"][k] = v

    for k, v in ABAQUS_IN["inputs"].items():
        if not "material" in k:  # only optimise on float type variables
            resp["input_data"]["design"][k] = v

    for k, v in ABAQUS_IN["outputs"].items():
        resp["output_data"]["design"][k] = v

    # copy input files to run folder
    analysis_input_folder_path = Path(resp["parameters"]["analysis_input_folder_path"])
    if not analysis_input_folder_path.is_dir():
        raise IsADirectoryError()
    for pattern in ["**/*.inp", "**/*.sim"]:
        for p in analysis_input_folder_path.glob(pattern):
            shutil.copy2(p, run_folder / p.name)

    # copy the post script (Note: from python 3.11.3 subprocess cannot call script in cwd)
    if "analysis_post_script" in resp["parameters"]:
        infile = Path(resp["parameters"]["analysis_post_script"])
    else:
        infile = Path("extract_constraints_data_from_odb.py")
    shutil.copy2(infile, run_folder / infile.name)

    # if there are RTS plies defined, read the shell mesh and save to setup data
    for rts_type in ["rts_plies", "rts_discrete_reinforcements"]:
        if rts_type in resp["parameters"]:
            for part in resp["parameters"][rts_type]:
                rts.setup(part, resp["parameters"], run_folder)

    # delete previous buckling cache file
    if (run_folder / "buckling_cache.json").is_file():
        os.remove(run_folder / "buckling_cache.json")

    resp["message"] = f"{datetime.now().strftime('%Y%m%d-%H%M%S')}: Setup completed."

    global SETUP_DATA
    SETUP_DATA = resp
    # print(SETUP_DATA)

    return resp


def compute(msg):
    resp = deepcopy(SETUP_DATA)
    run_folder = resp["parameters"]["outputs_folder_path"]
    resp.pop("input_data")

    ply_thickness_ratios = {}

    for part in resp["parameters"]["composite"]:
        # update all straight ply variables first
        # do this for all parts
        c_inputs = {
            k: v
            for k, v in msg["inputs"]["design"].items()
            if k.startswith(f"composite.{part}.") and not ".RTS." in k
        }
        for k, v in c_inputs.items():
            if isinstance(v, list):
                v = v[0]
            kchain = k.split(".")
            layup = kchain[2]
            ply = int(kchain[3])
            v_type = kchain[4]
            if v_type in ["thickness", "angle"]:
                resp["parameters"]["composite"][part][layup]["plies"][ply][v_type] = (
                    float(v)
                )
            else:
                raise ValueError(
                    f"Cannot optimise on non-float variable of type {v_type}."
                )

        def update_rts_parameters_with_design_inputs(rts_type):
            # first update resp["parameters"]["rts_plies"] with input values
            c_inputs = {
                k: v
                for k, v in msg["inputs"]["design"].items()
                if k.startswith(f"composite.{part}.RTS.")
            }
            for k, v in c_inputs.items():
                if isinstance(v, list):
                    v = v[0]
                kchain = k.split(".")
                ply_name = kchain[3]
                ctrl_point = int(kchain[4])
                v_type = kchain[5]
                if v_type in ["d", "inc_angle"]:
                    resp["parameters"][rts_type][part][ply_name]["control_pts"][
                        ctrl_point
                    ][v_type] = float(v)
                else:
                    raise ValueError(
                        f"Cannot optimise on non-float variable of type {v_type}."
                    )

        # TODO - I think this may fail if there are both complete rts plies and discrete reinforcements
        if (
            "rts_plies" in resp["parameters"]
            and part in resp["parameters"]["rts_plies"]
        ):
            update_rts_parameters_with_design_inputs("rts_plies")
            # write rts laminate inputs
            ply_thickness_ratios.update(
                rts.compute(part, resp["parameters"], run_folder)
            )
        elif (
            "rts_discrete_reinforcements" in resp["parameters"]
            and part in resp["parameters"]["rts_discrete_reinforcements"]
        ):
            update_rts_parameters_with_design_inputs("rts_discrete_reinforcements")
            # write rts laminate inputs
            ply_thickness_ratios.update(
                rts.compute(part, resp["parameters"], run_folder)
            )
        else:
            # write UD material parts input data to file
            ply_thickness_ratios.update(
                write_composite_inputs(
                    run_folder, part, resp["parameters"]["composite"][part]
                )
            )

    # delete previous output if it exists
    # if (run_folder / "odb_output.json").is_file():
    #     os.remove(run_folder / "odb_output.json")

    # execute abaqus
    while True:
        infile = run_folder / resp["parameters"]["analysis_file"]
        if (infile.parent / (infile.stem+".odb")).is_file():
            os.remove(infile.parent / (infile.stem+".odb"))
            print(f"deleted old odb @ {(infile.parent / (infile.stem+'.odb'))}\n")
        fea_resp = execute_fea(infile=infile, run_folder=run_folder)
        if not fea_resp["returncode"] == 0:
            raise ChildProcessError(
                f'abaqus returned non-zero exit status {resp["returncode"]}'
            )

        # read the ouputs from the odb
        if "analysis_post_script" in resp["parameters"]:
            infile = Path(resp["parameters"]["analysis_post_script"])
        else:
            infile = Path("extract_constraints_data_from_odb.py")

        post_resp = execute_fea_post(infile=infile, run_folder=run_folder)

        #### NOTE: custom output processing depending on output types
        # if "panel buckling modes" in post_resp["stdout"]:
        #     start_index = post_resp["stdout"].find("Only found") + 10
        #     end_index = post_resp["stdout"].find("panel buckling modes")
        #     _modes = [
        #         int(m)
        #         for m in post_resp["stdout"][start_index:end_index].strip().split("/")
        #     ]
        #     n_modes = (
        #         get_master_modes(run_folder / resp["parameters"]["analysis_file"])
        #         + (_modes[1] - _modes[0]) * 5
        #     )

        #     with open(run_folder / "abaqus_restart.log", "a") as f:
        #         f.write(
        #             f"Abaqus only found {_modes[0]} out of {_modes[1]} buckling modes @ {datetime.now().strftime('%Y%m%d-%H%M%S')}\nTrying to run Abaqus again with {n_modes} modes ... \n "
        #         )
        #     update_master_with_modes(
        #         n_modes, file=run_folder / resp["parameters"]["analysis_file"]
        #     )
        #     continue
        if not "GOT RESULTS FROM ODB" in post_resp["stdout"]:
            print(post_resp["stdout"])
            with open(run_folder / "abaqus_crash.log", "a") as f:
                f.write(
                    f"Abaqus post failed @ {datetime.now().strftime('%Y%m%d-%H%M%S')}\nTrying to run Abaqus again in 10s ... \n"
                )
            sleep(10)
            continue  # try running abaqus again
            # raise ChildProcessError(
            #     f'abaqus post returned non-zero exit status {post_resp["returncode"]}'
            # )
        else:
            break

    with open(run_folder / "odb_output.json", "r") as f:
        odb_data = json.load(f)

    # parse the abaqus outputs object into a local dictionary
    outputs = resp.pop("output_data")

    # map outputs to json dictionary items:
    for k, v in outputs["design"].items():
        kchain = k.split(".")
        output_name= kchain[0]
        if output_name == 'MASS': 
            step_name = kchain[1]
            instance_name = kchain[2]
            elset_name = kchain[3] 
            v = OUTPUT_type[output_name](odb_data[output_name][step_name][instance_name][elset_name])
            outputs["design"][k] = v    
        elif output_name == 'U':
            step_name = kchain[1]
            instance_name = kchain[2]
            nset_name = kchain[3] 
            dof_name = kchain[4]
            v = OUTPUT_type[output_name](odb_data[output_name][step_name][instance_name][nset_name][dof_name])
            outputs["design"][k] = v  
        elif output_name == 'T-RATIO-MIN' or output_name == 'T-RATIO-MAX': 
            instance_name = kchain[1] 
            elset_name = kchain[2].upper()
            ply_number = int(kchain[3])
            ratio_keys = {"T-RATIO-MIN": "min", "T-RATIO-MAX": "max"}  
            v = OUTPUT_type[output_name](ply_thickness_ratios[ratio_keys[output_name]][elset_name][ply_number])
            outputs["design"][k] = v
        elif output_name == 'TSAIW':
            step_name = kchain[1]
            instance_name = kchain[2]
            elset_name = kchain[3] 
            ply_number = kchain[4]
            v = odb_data[output_name][step_name][instance_name][elset_name][ply_number]
            outputs["design"][k] = OUTPUT_type[output_name](v["max"]) 
        elif output_name == 'E':
            step_name = kchain[1]
            instance_name = kchain[2]
            elset_name = kchain[3] 
            ply_number = kchain[4]
            strain_type = kchain[5]
            v = odb_data[output_name][step_name][instance_name][elset_name][ply_number][strain_type]
            outputs["design"][k] = OUTPUT_type[output_name](v["max"])   
       
            
    response = {}
    response["outputs"] = outputs
    response["message"] = (
        f"{datetime.now().strftime('%Y%m%d-%H%M%S')}: Compute completed."
    )

    print(response["message"])

    return response


def update_master_with_modes(n_modes, file):
    with open(file, "r") as f:
        data = f.readlines()

    is_buckle_keyword = False
    for ii, l in enumerate(data):
        if "*BUCKLE" in l.upper():
            is_buckle_keyword = True
            continue
        if is_buckle_keyword:
            options = l.split(",", 1)
            data[ii] = f"{n_modes},{options[1]}"
            break
    if not is_buckle_keyword:
        raise ValueError(f"Couldn't find the *Buckle keyword in input file to update.")

    with open(file, "w") as f:
        f.writelines(data)


def get_master_modes(file):
    with open(file, "r") as f:
        data = f.readlines()

    is_buckle_keyword = False
    for ii, l in enumerate(data):
        if "*BUCKLE" in l.upper():
            is_buckle_keyword = True
            continue
        if is_buckle_keyword:
            options = l.split(",", 1)
            break
    if not is_buckle_keyword:
        raise ValueError(f"Couldn't find the *Buckle keyword in input file to update.")

    return int(options[0])


def write_composite_inputs(run_folder, part, composites, ply_sectionPoints=3):
    sections = ""
    ply_thickness_ratios = {}
    for layup_name, layup_val in composites.items():
        elset = layup_val["elset"].upper()
        ori = layup_val["orientation"]
        offset = layup_val["offset"]
        sections += f"*Shell Section, elset={elset}, composite, orientation={ori}, offset={offset}, density=0., layup={layup_name}\n"

        t_layup = 0
        t_plies = [0] * len(layup_val["plies"])
        for ii, p in enumerate(layup_val["plies"]):
            sections += f"{p['thickness']:0.6f}, {ply_sectionPoints:d}, {p['material']}, {p['angle']:0.3f}, {p['name']}\n"
            if not p["name"] in PLIES_EXCLUDED_FROM_THICKNESS_RATIO:
                t_layup += p["thickness"]
                t_plies[ii] = p["thickness"]
        sections += "***\n"

        ply_thickness_ratios[elset] = [p / t_layup for p in t_plies]

    with open(run_folder / (part + "_sections.inp"), "w") as f:
        f.write(sections)

    return ply_thickness_ratios


def track_modes_by_MAC_criterion(
    buckling,
    run_folder: Path,
    cache="buckling_cache.json",
    retain_n_modes=None,
    threshold=0.7,  # MAC coefficients threshold below which no match is found
):
    """Track modes in buckling output by:
    1) performing a MAC analysis with the reference modes, and
    2) sorting the buckling modes to provide best correlation with previous modes (if available)
    """

    # check if this is a major iteration (with gradient calc) or not
    is_major_iter, is_line_search, is_finite_diff_step = is_major_iteration(run_folder)

    # cache the buckling data if this is the first analysis of the optimisation run or a major iteration
    if is_major_iter or is_line_search or not (run_folder / cache).is_file():
        buckling_sorted = {}
        for key in ["eigenvalues", "panel_ids", "modes"]:
            buckling_sorted[key] = {
                k: v
                for k, v in buckling[key].items()
                if int(k) in range(retain_n_modes)
            }
        set_cache_buckling(buckling_sorted, run_folder, cache)
    else:
        # read reference modes from cache
        ref_buckling = read_cache_buckling(run_folder, cache)
        ref_mode_order = list(ref_buckling["eigenvalues"].keys())

        # loop through every pair of modes and calculate the MAC criterion for every pair (optionally: plot)
        mac = get_mac(
            buckling,
            ref_buckling,
            plot_flag=PLOT_OPTION,
            figure_name=str(Path(run_folder, "mac.png")),
            threshold=threshold,
        )

        # track modes
        tracked_mode_order = sort_tracked_modes(
            mac,
            threshold=threshold,
        )
        if not tracked_mode_order == ref_mode_order:
            print(f"Found buckling mode switch: {tracked_mode_order}")

        buckling_sorted = {}
        for key in ["eigenvalues", "panel_ids", "modes"]:
            buckling_sorted[key] = {}
            for ii, mode in enumerate(ref_mode_order):
                mode_key = tracked_mode_order[ii]
                if mode_key == "not_tracked":
                    # ASSUMPTION NOTE: default to unchanged mode properties from the previous analysis
                    buckling_sorted[key][mode] = ref_buckling[key][mode]
                else:
                    buckling_sorted[key][mode] = buckling[key][mode_key]

    return list(buckling_sorted["eigenvalues"].values())


def set_cache_buckling(buckling, run_folder, cache):
    with open(run_folder / cache, "w") as f:
        json.dump(buckling, f)

    return None


def is_major_iteration(run_folder):
    """This is a workaround to determine the major iterations by reading the driver log file, if it exists.
    returns:  is_major_iter, is_line_search, is_finite_diff_step
    """

    if (run_folder / "run_driver.log").is_file():
        slsqp_iteration_start = False
        slsqp_iteration_completed = False

        with open(run_folder / "run_driver.log", "r") as f:
            data = f.readlines()

        data.reverse()
        for line in data:
            if (
                line.startswith("Driver total derivatives for iteration:")
                and slsqp_iteration_start
            ):
                # we are in a major iteration step
                return True, False, False
            elif (
                line.startswith("Driver total derivatives for iteration:")
                and not slsqp_iteration_start
            ) or (slsqp_iteration_completed and slsqp_iteration_start):
                # we are in a finite difference calculation step
                return False, False, True
            elif (
                line.startswith("Driver debug print for iter coord:")
                and not slsqp_iteration_start
            ):
                slsqp_iteration_start = True
            elif line.startswith("Objectives") and not slsqp_iteration_start:
                slsqp_iteration_completed = True
            elif (
                line.startswith("Driver debug print for iter coord:")
                and slsqp_iteration_start
            ):
                # 2 iteration starts in a row: we are in a line search step
                return False, True, False

        # This is the SLSQP|0 iteration
        return True, False, False
    else:
        # we are not running an optimisation
        return True, False, False


def read_cache_buckling(run_folder, cache):
    with open(run_folder / cache, "r") as f:
        ref_buckling = json.load(f)

    return ref_buckling


def get_mac(
    b,
    ref_b,
    plot_flag,
    figure_name,
    threshold,
):
    """Calculate the Modal assurance criterion (MAC) to identify mode switches."""

    def get_eigenvectors(modes):
        nodes = None
        eigv = []
        for mode in modes.values():
            n = [int(k) for k in mode.keys()]
            if not nodes:
                nodes = n
            elif not nodes == n:
                raise ValueError(
                    "Node ids in eigenvectors across modes of same model are not consistent."
                )  # not handled
            eigv.append(list(mode.values()))
        return np.array(eigv), np.array(nodes)

    def sort_by_nodes(mu_model2, nodes2, nodes1):
        # sort in case new and reference mode node numbering is not consistent
        if all(nodes1 == nodes2):
            return mu_model2
        else:
            mu_sorted = np.zeros_like(mu_model2)
            for ii, n in enumerate(nodes1):
                mu_sorted[:, ii] = mu_model2[:, np.where(nodes2 == n)[0][0]]
            return mu_sorted

    # reference modes
    n_modes_1 = len(ref_b["eigenvalues"])
    mu_model1, nodes1 = get_eigenvectors(ref_b["modes"])

    # new modes
    n_modes_2 = len(b["eigenvalues"])
    mu_model2, nodes2 = get_eigenvectors(b["modes"])
    mu_model2 = sort_by_nodes(mu_model2, nodes2, nodes1)

    # calculate the MACXP criterion from equation [28] in the reference
    mac = np.zeros([n_modes_1, n_modes_2])
    for ii, mii in enumerate(mu_model1):
        for jj, mjj in enumerate(mu_model2):
            coef = np.dot(mii, mjj) ** 2 / (np.dot(mii, mii) * np.dot(mjj, mjj))
            mac[ii, jj] = coef

    if plot_flag:
        ax = plt.figure().add_subplot()
        im = ax.imshow(mac, cmap="jet", vmin=0.0, vmax=1.0)
        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        colorbar = plt.colorbar(im, cax=cax)
        # Loop over data dimensions and create text annotations.
        for i in range(n_modes_1):
            for j in range(n_modes_2):
                if mac[i, j] >= threshold:
                    text = ax.text(
                        j,
                        i,
                        "%.2f" % mac[i, j],
                        ha="center",
                        va="center",
                        color="w",
                        fontsize=8,
                    )
        ax.set_title(f"MAC criterion")
        plt.xlabel("new modes")
        plt.ylabel("reference modes")
        if figure_name:
            plt.savefig(figure_name)
        # plt.show(block=True) # for debugging
        plt.close()

    return mac


def sort_tracked_modes(mac, threshold):
    # match eigs by highest mac coefficient
    tracked_modes_order = []
    for ii in range(mac.shape[0]):
        index = np.argmax(mac[ii, :])
        if mac[ii, index] >= threshold:
            tracked_modes_order.append(str(index))
        else:
            # if no match to ref can be found then set to NA
            tracked_modes_order.append("not_tracked")

    return tracked_modes_order


def execute_fea(infile: Path, run_folder: Path):
    """Run abaqus to generate the FEA output files."""

    if LOCAL_EXECUTES["ABAQUS"]:
        resp = subprocess.run(
            LOCAL_EXECUTES["ABAQUS"] + " job=" + str(infile.stem),
            cwd=run_folder,
            shell=True,
            check=False,
            capture_output=True,
        )
        return {"stdout": resp.stdout.decode("ascii"), "returncode": resp.returncode}
    else:
        raise ValueError("Need to specify an execution path for ABAQUS.")


def execute_fea_post(infile: Path, run_folder: Path):
    """Run abaqus to generate the FEA output files."""

    if LOCAL_EXECUTES["ABAQUS_POST"]:
        resp = subprocess.run(
            LOCAL_EXECUTES["ABAQUS_POST"] + " " + str(infile.name),
            cwd=run_folder,
            shell=True,
            check=False,
            capture_output=True,
        )
        return {"stdout": resp.stdout.decode("ascii"), "returncode": resp.returncode}
    else:
        raise ValueError("Need to specify an execution path for ABAQUS.")


if __name__ == "__main__":
    json_input_file = "test_dummy_master_linked_variables.json"

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

    resp = compute(TEST_MSG_COMPUTE)
    print(resp)
    with open(json_input_file[:-5] + "_outputs.log", "w") as f:
        f.write(resp["message"])

    print("TEST abaqus.py passed.")
