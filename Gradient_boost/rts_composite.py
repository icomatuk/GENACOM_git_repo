""" Utility functions for RTS composite modelling and analysis."""

from pathlib import Path
from math import sqrt, floor, ceil, cos, sin, radians, log10
import numpy as np
import copy
from shutil import copy2
import warnings
from scipy.spatial.transform import Rotation as R
import json

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import time

PLOT_OPTION = False
MODEL_CHECK_OPTION = False
TEST_PARAMETERS = {
    "rts_mesh_files": {
        "skin": {
            "mesh": "element_2_skin_baseline.inp",
            "shell_sections": "skin_sections.inp",
            "element_sets": "skin_element_sets.inp",
        }
    },
    "rts_plies": {
        "skin": {
            "SKIN_Ply-5": {
                "control_pts": [
                    {"d": -250.0, "id": "p1", "inc_angle": 45.0},
                    {"d": 0, "id": "p1", "inc_angle": 45.0},
                    {"d": 300.0, "id": "p2", "inc_angle": 0.0},
                    {"d": 1250.0, "id": "p2", "inc_angle": 0.0},
                ],
                "path_end": [-763.939, -824.999, 137.982],
                "path_start": [-1076.835, 0.000, 164.462],
                "orientations": {},
            }
        }
    },
    "rts_discrete_reinforcements": {
        "skin": {
            "SKIN_Ply-3": {
                "control_pts": [
                    {"d": -250.0, "id": "p1", "inc_angle": 0.0},
                    {"d": 0, "id": "p1", "inc_angle": -45.0},
                    {"d": 300.0, "id": "p2", "inc_angle": 0.0},
                    {"d": 1250.0, "id": "p2", "inc_angle": 45.0},
                ],
                "orientations": {},
                "path_end": [-763.939, -824.999, 137.982],
                "path_start": [-1076.835, 0.000, 164.462],
                "width": 50.0,
                "xinc": 10.0,
                "y_offset": 100.0,
            },
            "SKIN_Ply-4": {
                "control_pts": [
                    {"d": -60.0, "id": "p-start-fixed", "inc_angle": -10.0},
                    {"d": 40.0, "id": "p1", "inc_angle": -10.0},
                    {"d": 140.0, "id": "p2", "inc_angle": -10.0},
                    {"d": 240.0, "id": "p3", "inc_angle": -10.0},
                    {"d": 340.0, "id": "p4", "inc_angle": -10.0},
                    {"d": 440.0, "id": "p5", "inc_angle": -10.0},
                    {"d": 660.0, "id": "p-end-fixed", "inc_angle": -10.0},
                ],
                "path_end": [-1342.723713, -591.170688, 112.958736],
                "path_start": [-1076.835, 0.000, 164.462],
                "orientations": {},
                "width": 20.0,
                "xinc": 10.0,
                "y_offset": -20.0,
            },
        }
    },
    "composite": {
        "skin": {
            "ELE2_SKIN_TOP_0": {
                "elset": "top_0",
                "orientation": "CSYS-1-SKIN",
                "offset": "SNEG",
                "plies": [
                    {
                        "angle": 45.0,
                        "thickness": 0.209,
                        "material": "PC_IT-E49-3_BIAX",
                        "name": "SKIN_Ply-1",
                    },
                    {
                        "angle": 45.0,
                        "thickness": 0.205,
                        "material": "PC_HT-E49-3_BIAX",
                        "name": "SKIN_Ply-2",
                    },
                    {
                        "angle": 0.0,
                        "thickness": 0.205,
                        "material": "PC_IU-E49-3_UD",
                        "name": "SKIN_Ply-3",
                    },
                    {
                        "angle": 45.0,
                        "thickness": 0.205,
                        "material": "PC_IU-E49-3_UD",
                        "name": "SKIN_Ply-4",
                    },
                    {
                        "angle": -45.0,
                        "thickness": 0.205,
                        "material": "PC_IU-E49-3_UD",
                        "name": "SKIN_Ply-5",
                    },
                ],
            },
            "ELE2_SKIN_TOP_1": {
                "elset": "top_1",
                "orientation": "CSYS-1-SKIN",
                "offset": "SNEG",
                "plies": [
                    {
                        "angle": 45.0,
                        "thickness": 0.209,
                        "material": "PC_IT-E49-3_BIAX",
                        "name": "SKIN_Ply-1",
                    },
                    {
                        "angle": 45.0,
                        "thickness": 0.205,
                        "material": "PC_HT-E49-3_BIAX",
                        "name": "SKIN_Ply-2",
                    },
                    {
                        "angle": 0.0,
                        "thickness": 0.205,
                        "material": "PC_IU-E49-3_UD",
                        "name": "SKIN_Ply-3",
                    },
                    {
                        "angle": 45.0,
                        "thickness": 0.205,
                        "material": "PC_IU-E49-3_UD",
                        "name": "SKIN_Ply-4",
                    },
                    {
                        "angle": -45.0,
                        "thickness": 0.205,
                        "material": "PC_IU-E49-3_UD",
                        "name": "SKIN_Ply-5",
                    },
                ],
            },
        }
    },
}

IGNORE_ELSETS = [
    # any non-design elsets, for example sets of sets that are only used for post-processing
]

ANGLE_ACCURACY = 1
PLY_NB_ACCURACY = 10
PLIES_EXCLUDED_FROM_THICKNESS_RATIO = [
    # e.g. outer woven pies that are not being designed and where thickness ratio constraints do not apply
]

SETUP_DATA = {}


def setup(part, inputs, run_folder):
    elements, nodes, baseline_elsets, coordinate_systems, flat_coordinate_systems = (
        get_mesh(inputs["rts_mesh_files"][part], run_folder)
    )

    ori_dist = {}
    if inputs["rts_mesh_files"][part].get("flattened_mesh"):

        for ori_name, coord in coordinate_systems.items():

            ori_dist_name = "dist_" + ori_name
            if not ori_dist.get(ori_dist_name):
                ori_dist[ori_dist_name] = {}

            for eid, ele in elements.items():
                ori_dist[ori_dist_name][eid] = (
                    ele["flattened"]["r_flat_to_3D"].as_matrix() @ coord["R"].T
                )

        with Timer("write_ori_distributions"):
            write_orientation_ditributions(
                run_folder,
                file="orientation_distributions.inp",
                ori_dist=ori_dist,
                default_ori=coordinate_systems,
            )
        write_distribution_tables(
            run_folder,
            file="distribution_tables.inp",
            ori_dist=ori_dist,
        )

    global SETUP_DATA
    SETUP_DATA[part] = {
        "nodes": nodes,
        "elements": elements,
        "baseline_elsets": baseline_elsets,
        "coordinate_systems": coordinate_systems,
        "flat_coordinate_systems": flat_coordinate_systems,
        "ori_dist": ori_dist,
    }

    print(f"Completed rts setup for part {part}.")

    return None


def compute(part, inputs, run_folder):
    # get elements and baseline sets from setup data
    elements = SETUP_DATA[part]["elements"]
    baseline_elsets = SETUP_DATA[part]["baseline_elsets"]

    # get rts ply property functions for all plies
    if "rts_plies" in inputs:
        rts_plies = inputs["rts_plies"][part]
        for ply in rts_plies:
            set_rts_distribution(rts_plies[ply], ANGLE_ACCURACY)
    else:
        rts_plies = {}

    if "rts_discrete_reinforcements" in inputs:
        rts_discrete_reinforcements = inputs["rts_discrete_reinforcements"][part]
        for ply in rts_discrete_reinforcements:
            set_rts_distribution(
                rts_discrete_reinforcements[ply],
                ANGLE_ACCURACY,
                rts_discrete_reinforcements=True,
                name=ply,
                run_folder=run_folder,
            )
    else:
        rts_discrete_reinforcements = {}

    # lookup new properties for all elements
    element_sets, shell_sections, ply_thickness_ratios = get_element_shell_sections(
        laminates=inputs["composite"][part],
        baseline_elsets=baseline_elsets,
        elements=elements,
        rts_plies=rts_plies,
        rts_discrete_reinforcements=rts_discrete_reinforcements,
        angle_accuracy=ANGLE_ACCURACY,
        coordinate_systems=SETUP_DATA[part]["coordinate_systems"],
    )

    # check element vs applied CS orientations
    if MODEL_CHECK_OPTION == True:
        ### compare the projection of the steering path to the material reference direction - these should be identical within a small tol
        if rts_plies:
            check_rts_path_angles(rts_plies, inputs["composite"][part])
        if rts_discrete_reinforcements:
            check_rts_path_angles(
                rts_discrete_reinforcements, inputs["composite"][part]
            )
        # check element normals wrt. local coordinate system z-direction
        check_shell_element_normals(
            element_sets,
            elements,
            shell_sections,
            nodes=SETUP_DATA[part]["nodes"],
            coordinate_systems=SETUP_DATA[part]["coordinate_systems"],
        )

    # plot element centroids by property
    if PLOT_OPTION == True:
        view_projection = None
        if inputs == TEST_PARAMETERS:
            view_projection = {"elev": -90, "azim": 0, "roll": 180}
        plot_3D_scatter_properties(
            element_sets,
            elements,
            shell_sections,
            rts_discrete_reinforcements=rts_discrete_reinforcements,
            view_projection=view_projection,
        )

    if inputs["rts_mesh_files"][part].get("flattened_mesh"):
        # for each shell section orientation, replace orientation with a distribution
        # the distribution default is the original orientation definition
        # then for all elements the value is the original orientation definition rotated
        # by r in the element flattened dictionary - moved to SETUP

        for _, layup_val in shell_sections.items():
            ref_ori_name = layup_val["orientation"].upper()
            ori_dist_name = "dist_" + ref_ori_name
            layup_val["orientation"] = ori_dist_name

    # write elsets and shell sections to file
    write_element_sets(
        element_sets,
        run_folder=run_folder,
        file=inputs["rts_mesh_files"][part]["element_sets"],
    )
    write_composite_inputs(
        run_folder,
        shell_sections,
        file=inputs["rts_mesh_files"][part]["shell_sections"],
        elsets_file=inputs["rts_mesh_files"][part]["element_sets"],
        orientation_ditributions="orientation_distributions.inp",
        angle_accuracy=ANGLE_ACCURACY,
        ori_dist=SETUP_DATA[part]["ori_dist"],
    )

    print(
        f"Completed writing of RTS section input data for part {part}.\n Number of composite shell sections in {part} = {len(shell_sections)}."
    )

    return ply_thickness_ratios


def get_mesh(mesh_files, run_folder):
    # read the 2D mesh from CGX output file
    elements, nodes, baseline_elsets, coordinate_systems, flat_coordinate_systems = (
        get_2Dmesh(mesh_files, run_folder)
    )

    # 2) Calculate the element normals, and then the node normal
    corner_nodes = set()
    #### set_node_normals(nodes=nodes, elements=elements)
    set_element_centroid(elements=elements, nodes=nodes, corner_nodes=corner_nodes)

    # get local coordinate_systems rotation matrices for fast vector transforms
    if coordinate_systems:
        get_transformation_matrices(coordinate_systems)

    if flat_coordinate_systems:
        get_transformation_matrices(flat_coordinate_systems)

    return elements, nodes, baseline_elsets, coordinate_systems, flat_coordinate_systems


def get_2Dmesh(mesh_files, run_folder):
    elements, nodes, baseline_elsets, coordinate_systems = get_nodes_and_elements(
        run_folder, file=mesh_files["mesh"]
    )

    flat_coordinate_systems = None
    if "flattened_mesh" in mesh_files:

        print("reading flattened mesh from file")
        flat_elements, flat_nodes, _, flat_coordinate_systems = get_nodes_and_elements(
            run_folder, file=mesh_files["flattened_mesh"]["file"]
        )

        # make sure the reference element location is a match in both meshes
        # otherwise rotate and translate flattened mesh + plot overlay
        check_flattened_mesh(
            flat_nodes=flat_nodes,
            nodes=nodes,
            elements=elements,
            flat_elements=flat_elements,
            reference_element_id=mesh_files["flattened_mesh"]["reference_element_id"],
        )

        # for each element calculate the rotation from the flat to the curved shape
        # and save rotation matrix and flattened element centroid xyz to the elements dictionary
        set_element_centroid(
            elements=flat_elements, nodes=flat_nodes, corner_nodes=set()
        )

        for id, ele in elements.items():
            rotation, _, _ = get_transforms(
                elements[id],
                nodes,
                flat_nodes,
                flat_elements[id],
            )
            ele["flattened"] = {
                "r_flat_to_3D": rotation,
                "flat_xyz": flat_elements[id]["global_xyz"],
            }

    return (
        elements,
        nodes,
        baseline_elsets,
        coordinate_systems,
        flat_coordinate_systems,
    )


def get_nodes_and_elements(run_folder, file):
    with open(run_folder / file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    read_nodes = False
    read_elements = False
    read_elset = False
    read_coord = False
    ele_type = None
    nodes = {}
    elements = {}
    baseline_elsets = {}
    coordinate_systems = {}
    for line in lines:
        # --------------------------------------------------------------------------------
        # set read flags and parameters if a new keyword is detected
        if line.strip().upper().startswith("*NODE"):
            print("start reading nodes.")
            read_nodes = True
            read_elements = False
            read_coord = False
            read_elset = False
            ele_type = None
            continue
        elif line.strip().upper().startswith("*ELEMENT"):
            print("start reading elements.")
            read_nodes = False
            read_elements = True
            read_coord = False
            read_elset = False
            ele_params = {
                data.split("=")[0].strip(): data.split("=")[1].strip()
                for data in line.split(",")
                if "=" in data
            }
            ele_type = [{k: v} for k, v in ele_params.items() if k.upper() == "TYPE"][0]
            continue
        elif line.strip().upper().startswith("*ORIENTATION"):
            print("start reading orientation.")
            read_nodes = False
            read_elements = False
            read_coord = True
            read_elset = False
            ele_type = None

            coord_params = {
                data.split("=")[0].strip().upper(): data.split("=")[1].strip().upper()
                for data in line.split(",")
                if "=" in data
            }
            coord_name = coord_params["NAME"]
            if (
                "DEFINITION" in coord_params
                and not coord_params["DEFINITION"] == "COORDINATES"
            ):
                raise ValueError(
                    "Invalid option for DEFINITION parameter on *Orientation card found. Only DEFINITION=COORDINATES currently implemented."
                )
            if "SYSTEM" in coord_params and not coord_params["SYSTEM"] == "RECTANGULAR":
                raise ValueError(
                    "Invalid option for SYSTEM parameter on *Orientation card found. Only SYSTEM=RECTANGULAR currently implemented."
                )
            coordinate_systems[coord_name] = {}
            continue
        elif line.strip().upper().startswith("*ELSET"):
            read_nodes = False
            read_elements = False
            read_coord = False
            elset_params = {
                data.split("=")[0].strip().upper(): data.split("=")[1].strip().upper()
                for data in line.split(",")
                if "=" in data
            }
            elset_name = elset_params["ELSET"]
            if not elset_name in IGNORE_ELSETS:
                print(f"start reading element set {elset_name}.")
                read_elset = True
                if "GENERATE" in line.upper():
                    elset_generate = True
                else:
                    elset_generate = False
            else:
                print(f"this element set is ignored: {elset_name}")
            continue
        elif line.strip().startswith("**"):
            continue
        elif line.strip().startswith("*"):
            read_nodes = False
            read_elements = False
            read_coord = False
            read_elset = False
            continue

        # ---------------------------------------------------------------------------------
        # read data

        if read_nodes:
            data = line.strip().split(",")
            # parse nodes into { ID: {"global_xyz": [x,y,z]}}
            nodes[int(data[0])] = {"global_xyz": [float(v) for v in data[1:]]}
            continue

        if read_elements:
            data = line.strip().split(",")
            # parse elements into { ID: {"type": S8, "nodes": [1-8]}}
            elements[int(data[0])] = {
                "type": ele_type,
                "nodes": [int(v) for v in data[1:]],
            }
            continue

        if read_coord:
            data = line.strip().split(",")
            if not coordinate_systems[coord_name]:
                # first line
                if len(data) == 6:
                    data += [
                        "0.0",
                        "0.0",
                        "0.0",
                    ]  # The default location of the origin, c, is the global origin.
                pts = [float(v) for v in data]
                coordinate_systems[coord_name] = {
                    "a": pts[:3],
                    "b": pts[3:6],
                    "c": pts[6:],
                }
            else:
                # second line
                coordinate_systems[coord_name]["rotation_dir"] = int(data[0])
                try:
                    # angle
                    coordinate_systems[coord_name]["rotation_val"] = float(data[1])
                except:
                    # distribution name
                    coordinate_systems[coord_name]["rotation_val"] = data[1]

        if read_elset and elset_generate:
            data = line.strip().split(",")
            baseline_elsets[elset_name] = list(
                range(int(data[0]), int(data[1]) + 1, int(data[2]))
            )
            continue

        if read_elset and not elset_generate:
            data = line.strip().split(",")
            if elset_name in baseline_elsets:
                baseline_elsets[elset_name] += [int(d) for d in data if d]
            else:
                baseline_elsets[elset_name] = [int(d) for d in data if d]
            continue

    return elements, nodes, baseline_elsets, coordinate_systems


def check_flattened_mesh(
    flat_nodes, nodes, elements, flat_elements, reference_element_id, tol=1e-6
):
    r, t, r_centre = get_transforms(
        elements[reference_element_id],
        nodes,
        flat_nodes,
        flat_elements[reference_element_id],
    )

    if np.allclose(
        r.as_euler("xyz", degrees=True), np.array([0.0, 0.0, 0.0]), atol=tol
    ) and np.allclose(t, np.array([0.0, 0.0, 0.0]), atol=tol):
        pass
    else:
        print(
            f"transforming flattened mesh by R={r.as_matrix()} and T={t} to align reference element {reference_element_id}"
        )
        flat_nodes = transform_flatened(flat_nodes, r, t, r_centre)

    if PLOT_OPTION:
        # plot overlay of the flattened meshes
        plot_flat(nodes, elements, flat_nodes, flat_elements, target="global_xyz")

    return None


def transform_flatened(nodes, r, t, r_centre, transform_target="global_xyz"):
    transformed = {}
    for nid, val in nodes.items():
        xyz = val[transform_target]
        transformed[nid] = {
            "global_xyz": r.apply(np.array(xyz) - r_centre) + r_centre + t
        }
    return transformed


def get_transforms(ele, nodes, flat_nodes, flat_ele):
    nids = ele["nodes"]
    nodes_ref = [nodes[nid]["global_xyz"] for nid in nids]
    ab_ref = points2unit_vector(nodes_ref[1], nodes_ref[0])
    ac_ref = points2unit_vector(nodes_ref[2], nodes_ref[0])

    r_centre = [0.0, 0.0, 0.0]
    for xyz in nodes_ref:
        r_centre = [c1 + c2 for c1, c2 in zip(r_centre, xyz)]
    r_centre = np.array([c / len(ele["nodes"]) for c in r_centre])

    flat_nids = flat_ele["nodes"]
    nodes_flat = [flat_nodes[nid]["global_xyz"] for nid in flat_nids]
    ab_flat = points2unit_vector(nodes_flat[1], nodes_flat[0])
    ac_flat = points2unit_vector(nodes_flat[2], nodes_flat[0])

    # rotation from flat to reference coordinate system
    r = R.align_vectors(
        [
            ab_ref,
            ac_ref,
        ],
        [ab_flat, ac_flat],
    )[0]

    # translation from flat to reference coordinate system
    t = nodes_ref[0] - r.apply(nodes_flat[0])

    return r, t, r_centre


def plot_flat(nodes, elements, flat_nodes, flat_elements, target="global_xyz"):

    ref_nodes = np.empty([len(nodes), 3])
    f_nodes = np.empty([len(flat_nodes), 3])
    counter = 0
    ref_node_lookup = {}
    for nid, v in nodes.items():
        ref_node_lookup[nid] = counter
        ref_nodes[counter, :] = v[target]
        counter += 1

    f_nodes_lookup = {}
    counter = 0
    for nid, v in flat_nodes.items():
        f_nodes_lookup[nid] = counter
        f_nodes[counter, :] = v[target]
        counter += 1

    # ref_surface = [
    #     [ref_nodes[ref_node_lookup[nid]] for nid in e["nodes"]]
    #     for _, e in elements.items()
    # ]
    # flat_surface = [
    #     [f_nodes[f_nodes_lookup[nid]] for nid in e["nodes"]]
    #     for _, e in flat_elements.items()
    # ]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.view_init(elev=0, azim=0, roll=0)
    ax.scatter(ref_nodes[:, 0], ref_nodes[:, 1], ref_nodes[:, 2], marker="o")
    ax.scatter(f_nodes[:, 0], f_nodes[:, 1], f_nodes[:, 2], marker="x")
    # ax.add_collection3d(Poly3DCollection(ref_surface, color="b", edgecolor="k"))
    # ax.add_collection3d(Poly3DCollection(flat_surface, color="r", edgecolor="k"))

    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_zlabel("Z-axis")
    ax.set_aspect("equal", adjustable="box")
    # ax.legend(["actual mesh", "flattened mesh"])

    plt.show()


def set_element_centroid(elements, nodes, corner_nodes):
    for _, element in elements.items():
        centroid = [0.0, 0.0, 0.0]
        for nid in element["nodes"]:
            xyz = nodes[nid]["global_xyz"]
            centroid = [c1 + c2 for c1, c2 in zip(centroid, xyz)]
        centroid = [c / len(element["nodes"]) for c in centroid]
        element["global_xyz"] = centroid

        ## TODO for each element normal is the average of the normals calculated at the elemnent nodes from the element edges ( sum of vectors / nb nodes)
        # (assume n = [0,0,1] for the flat panel)
        # element["normal"] = [0.0, 0.0, 1.0]
        if element["type"] == "S8":
            corner_nodes.update(element["nodes"][:4])

    return None


def points2unit_vector(a, b):
    upath = np.array(a) - np.array(b)
    return upath / np.linalg.norm(upath)


def get_transformation_matrices(coordinate_systems):
    for name, coord in coordinate_systems.items():
        x = points2unit_vector(coord["a"], coord["c"])
        z = np.cross(x, points2unit_vector(coord["b"], coord["c"]))
        y = np.cross(z, x)
        # ["R"].T is the rotation from the local to the global CS
        # ["R"] is the rotation from the global to the local CS
        coordinate_systems[name]["R"] = np.array([x, y, z])

    return None


def set_node_normals(nodes, elements):
    for nid in nodes:
        ## TODO for each node, find the connected elements and average normal at the node as sum of vectors / nb connected elements
        # (assume n = [0,0,1] for the flat panel)
        nodes[nid]["normal"] = [0.0, 0.0, 1.0]

    return None


def set_rts_distribution(
    ply, angle_accuracy, rts_discrete_reinforcements=False, name=None, run_folder=None
):
    # get path unit vector
    upath = points2unit_vector(ply["path_end"], ply["path_start"])

    # linear interpolation
    ranges = []
    ctrl_points = []
    for p0, p1 in zip(ply["control_pts"][:-1], ply["control_pts"][1:]):
        ranges.append((p0["d"], p1["d"]))
        ctrl_points.append((p0["inc_angle"], p1["inc_angle"]))

    # input to f_inc_angle is dot product of upath with vector from path start to node
    ply["f_inc_angle"] = {
        "f": piecewise_linear_angle_interp,
        "args": [ranges, ctrl_points, angle_accuracy],
    }
    ply["upath"] = upath

    if rts_discrete_reinforcements:
        ply["tape_edges"] = get_tape_edge_paths(
            ply, plot_title=name, run_folder=run_folder
        )

    return None


def f_thickness_ratio(x):
    return 1 / cos(radians(x))


def piecewise_linear_angle_interp(x, ranges, ctrl_points, angle_accuracy=1):
    f_angles = lambda x, a0, a1, d0, d1: round(
        (a1 - a0) / (d1 - d0) * (x - d0) + a0, angle_accuracy
    )  #
    for d, a in zip(ranges, ctrl_points):
        if x >= d[0] and x < d[1]:
            return f_angles(x, a[0], a[1], d[0], d[1])
    raise ValueError(f"x value {x} is outside control point range.")


def get_tape_edge_paths(ply, plot_title, run_folder, plot=PLOT_OPTION):
    """calculate reinforcement tape boundary paths (ymin, ymax = ymin + w_tape)"""

    # integrate the fibre angle variations to obtain y=f(x) path definition
    # assuming linear fibre angle variations
    x = np.arange(ply["control_pts"][0]["d"], ply["control_pts"][-1]["d"], ply["xinc"])
    angles = np.array(
        [ply["f_inc_angle"]["f"](xii, *ply["f_inc_angle"]["args"]) for xii in x]
    )
    integration_angles = (angles[1:] - angles[:-1]) / 2 + angles[:-1]
    y = (
        np.cumsum(
            np.insert(np.tan(np.radians(integration_angles)) * (x[1:] - x[:-1]), 0, 0.0)
        )
        + ply["y_offset"]
    )

    tape_edges = [(x, y), (x, y + ply["width"])]

    if plot:
        ctrl_points = [p["d"] for p in ply["control_pts"]]
        ctrl_angles = [p["inc_angle"] for p in ply["control_pts"]]

        plt.plot(tape_edges[0][0], tape_edges[0][1], label="lower edge")
        plt.plot(tape_edges[1][0], tape_edges[1][1], label="upper edge")
        plt.plot(
            ctrl_points, np.zeros([len(ctrl_points)]), "ko", label="control points"
        )
        ax = plt.gca()
        ax.set_aspect("equal", adjustable="box")
        plt.title(plot_title)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.xlim([x[0], x[-1]])
        plt.grid(True)
        plt.legend(loc="right")
        plt.savefig(run_folder / f"rts_reinforcement_ply_{plot_title}.png")
        plt.show(block=True)

    return tape_edges


def update_laminate(
    laminate, ply_nb, material, xyz, rts_ply, set_name_suffix, angle_accuracy
):
    inc_angle = get_inc_angle_at_point(
        xyz,
        rts_ply["path_start"],
        rts_ply["upath"],
        rts_ply["f_inc_angle"],
    )
    thickness = f_thickness_ratio(inc_angle) * material["thickness"]
    total_angle = (
        inc_angle
        + rts_ply["orientations"][laminate["orientation"].upper()]["z_rot_angle"]
    )  # replaced material angle with actual calculated path 2D angle
    laminate["plies"][ply_nb]["angle"] = total_angle
    laminate["plies"][ply_nb]["thickness"] = thickness
    suffix = "_{angle:.{a_format}f}_{ply_nb:.{t_format}f}".format(
        angle=total_angle,
        a_format=angle_accuracy,
        ply_nb=ply_nb,
        t_format=PLY_NB_ACCURACY,
    )
    set_name_suffix += suffix.replace(".", "pt")

    return set_name_suffix


def get_element_shell_sections(
    laminates,
    baseline_elsets,
    elements,
    rts_plies,
    rts_discrete_reinforcements,
    angle_accuracy,
    coordinate_systems,
):
    def remove_plies_by_name(plies_to_remove, new_laminate):
        # remove the ply from the laminate definition by ply name
        while plies_to_remove:
            for jj, p in enumerate(new_laminate["plies"]):
                if p["name"] == plies_to_remove[0]:
                    new_laminate["plies"].pop(jj)
                    plies_to_remove.pop(0)
                    break
            continue

    element_sets = {}  # {elset_name: [eids], ...}
    shell_sections = {}  # same format as UD composite parameters
    set_ids = {}
    id_counter = 0
    ply_thickness_ratios = {"max": {}, "min": {}}
    for set_name, laminate in laminates.items():
        elset = baseline_elsets[laminate["elset"].upper()]
        ply_thickness_ratios["max"][laminate["elset"].upper()] = [0.0] * len(
            laminate["plies"]
        )
        ply_thickness_ratios["min"][laminate["elset"].upper()] = [0.0] * len(
            laminate["plies"]
        )
        ori = laminate["orientation"].upper()
        for eid in elset:
            element = elements[eid]
            if "flattened" in element:
                xyz = element["flattened"]["flat_xyz"]
            else:
                xyz = element["global_xyz"]
            new_laminate = copy.deepcopy(laminate)
            set_name_suffix = ""
            t_layup = 0
            t_plies = [0] * len(laminate["plies"])
            plies_to_remove = []
            for ii, ply in enumerate(laminate["plies"]):
                if ply["name"] in rts_plies:
                    if not ori in rts_plies[ply["name"]]["orientations"]:
                        get_ply_material_rotation_matrix(
                            rts_ply=rts_plies[ply["name"]],
                            ori=coordinate_systems[ori],
                            name=ori,
                        )
                    set_name_suffix = update_laminate(
                        new_laminate,
                        ply_nb=ii,
                        material=ply,
                        xyz=xyz,
                        rts_ply=rts_plies[ply["name"]],
                        set_name_suffix=set_name_suffix,
                        angle_accuracy=angle_accuracy,
                    )

                if ply["name"] in rts_discrete_reinforcements:
                    # determine if the element is within reinforcement or not
                    if (
                        not ori
                        in rts_discrete_reinforcements[ply["name"]]["orientations"]
                    ):
                        get_ply_material_rotation_matrix(
                            rts_ply=rts_discrete_reinforcements[ply["name"]],
                            ori=coordinate_systems[ori],
                            name=ori,
                        )
                    if is_element_in_reinforcement(
                        xyz,
                        rts_discrete_reinforcements[ply["name"]],
                        orientation=coordinate_systems[ori],
                        ori_name=ori,
                    ):
                        # replicate same functions as above
                        set_name_suffix = update_laminate(
                            new_laminate,
                            ply_nb=ii,
                            material=ply,
                            xyz=xyz,
                            rts_ply=rts_discrete_reinforcements[ply["name"]],
                            set_name_suffix=set_name_suffix,
                            angle_accuracy=angle_accuracy,
                        )
                    else:
                        plies_to_remove.append(new_laminate["plies"][ii]["name"])
                if (
                    not new_laminate["plies"][ii]["name"]
                    in PLIES_EXCLUDED_FROM_THICKNESS_RATIO
                    and not new_laminate["plies"][ii]["name"] in plies_to_remove
                ):
                    t_layup += new_laminate["plies"][ii]["thickness"]
                    t_plies[ii] = new_laminate["plies"][ii]["thickness"]
                elif (
                    not new_laminate["plies"][ii]["name"]
                    in PLIES_EXCLUDED_FROM_THICKNESS_RATIO
                    and new_laminate["plies"][ii]["name"] in plies_to_remove
                ):
                    t_layup += 0.0
                    t_plies[ii] = 0.0

            # TODO: it is possible to define the shell sections with distributions both for the fibre angle and for the
            # thickness (see work by Kostas) - use this approach in the future to reduce the number of shell sections!

            # Create sets of elements by ply and common inc_angle.
            elset_name = laminate["elset"].upper() + set_name_suffix
            elset_exists = elset_name in baseline_elsets
            if elset_name in set_ids:
                # add element to existing rts laminate element set
                element_sets[set_ids[elset_name]].append(eid)
            elif not elset_exists:
                # create new rts laminate shell section and element set
                id_counter += 1
                set_ids[elset_name] = laminate["elset"].upper() + f"_{id_counter}"
                section_name = set_name.upper() + f"_{id_counter}"

                remove_plies_by_name(plies_to_remove, new_laminate)
                new_laminate["elset"] = set_ids[elset_name]
                shell_sections[section_name] = new_laminate
                element_sets[set_ids[elset_name]] = [eid]
            else:
                # also generates shell sections for non-RTS laminates on the same part
                set_ids[elset_name] = laminate["elset"].upper() + "_0"
                section_name = set_name.upper() + "_0"

                remove_plies_by_name(plies_to_remove, new_laminate)
                new_laminate["elset"] = set_ids[elset_name]
                shell_sections[section_name] = new_laminate
                element_sets[set_ids[elset_name]] = [eid]

            if t_layup > 0.0:
                ply_thickness_ratios["max"][laminate["elset"].upper()] = [
                    p / t_layup if p / t_layup > ref else ref
                    for p, ref in zip(
                        t_plies, ply_thickness_ratios["max"][laminate["elset"].upper()]
                    )
                ]
                ply_thickness_ratios["min"][laminate["elset"].upper()] = [
                    p / t_layup if p / t_layup < ref or ref == 0.0 else ref
                    for p, ref in zip(
                        t_plies, ply_thickness_ratios["min"][laminate["elset"].upper()]
                    )
                ]
            else:
                raise ValueError(
                    f"Laminate {laminate['elset'].upper()} has zero total thickness!"
                )

    return element_sets, shell_sections, ply_thickness_ratios


def interpolate_2pts(x, y, x_val):
    return (y[1] - y[0]) / (x[1] - x[0]) * (x_val - x[0]) + y[0]


def get_ply_material_rotation_matrix(rts_ply, ori, name):
    # calculate rotation angle about local z-direction
    # ori["R"].T is the rotation from the local to the global CS
    # ori["R"] is the rotation from the global to the local CS
    upath_2d = (ori["R"] @ rts_ply["upath"])[:2]
    angle_2d = np.sign(upath_2d[1]) * np.rad2deg(
        np.arccos(upath_2d[0] / np.linalg.norm(upath_2d))
    )
    rts_ply["orientations"][name] = {"z_rot_angle": angle_2d}

    # calculate rotion matrix due to material rotation about z-direction
    matrix = R.from_euler("z", angle_2d, degrees=True).as_matrix()
    rts_ply["orientations"][name]["R"] = matrix


def is_element_in_reinforcement(xyz, ply, orientation, ori_name):
    # resolve vector in local coordinate system, only retain x, y to check against conditions
    material_R = ply["orientations"][ori_name]["R"]
    vector2node = (
        material_R.transpose()
        @ orientation["R"]
        @ (np.array(xyz) - np.array(ply["path_start"]))
    )

    edge_1 = ply["tape_edges"][0]
    edge_2 = ply["tape_edges"][1]
    if vector2node[0] >= edge_1[0][0] and vector2node[0] <= edge_1[0][-1]:
        # True if element is between control points (in x-direction)
        diff_x = np.absolute(edge_1[0] - vector2node[0])
        index = diff_x.argmin()
        if edge_1[0][index] - vector2node[0] > 0:
            ii = range(index - 1, index + 1)
        else:
            ii = range(index, index + 2)
        y1 = interpolate_2pts(x=edge_1[0][ii], y=edge_1[1][ii], x_val=vector2node[0])
        y2 = interpolate_2pts(x=edge_2[0][ii], y=edge_2[1][ii], x_val=vector2node[0])

        if vector2node[1] >= y1 and vector2node[1] <= y2:
            # True if element is also between reinforcement edges (in y-direction)
            return True

    return False


def get_inc_angle_at_point(xyz, start, upath, f_inc_angle):
    vector2node = np.array(xyz) - start
    projected_length = np.dot(upath, vector2node)
    inc_angle = f_inc_angle["f"](projected_length, *f_inc_angle["args"])
    return inc_angle


def write_element_sets(element_sets, run_folder=None, file="skin_element_sets.inp"):
    lines = []

    for name, elset in element_sets.items():
        lines.append(f"*ELSET,ELSET={name}\n")
        for chunk in divide_chunks(elset, 16):
            # note trailing comas are no a problem for ccx or abaqus
            lines.append(", ".join([f"{int(eid):d}" for eid in chunk]) + ",\n")

    # write string of input lines to file
    with open(run_folder / file, "w", encoding="utf-8") as f:
        f.write("".join(lines))


def divide_chunks(l, n):
    if not isinstance(l, list):
        l = list(l)
    # looping till length l
    for i in range(0, len(l), n):
        yield l[i : i + n]


def write_composite_inputs(
    run_folder,
    shell_sections,
    ply_sectionPoints=3,
    file="skin_sections.inp",
    elsets_file="skin_element_sets.inp",
    orientation_ditributions="orientation_distributions.inp",
    angle_accuracy=1,
    ori_dist=None,
):
    sections = f'*include, INPUT="{elsets_file}"\n***\n'

    if ori_dist:
        sections += f'*include, INPUT="{orientation_ditributions}"\n***\n'

    for layup_name, layup_val in shell_sections.items():
        elset = layup_val["elset"]
        ori = layup_val["orientation"]
        offset = layup_val["offset"]
        sections += f"*Shell Section, elset={elset}, composite, orientation={ori}, offset={offset}, density=0., layup={layup_name}, symmetric\n"

        for p in layup_val["plies"]:
            sections += f"{p['thickness']:0.6f}, {ply_sectionPoints:d}, {p['material']}, {p['angle']:0.{angle_accuracy}f}, {p['name']}\n"
        sections += "***\n"

    with open(run_folder / file, "w", encoding="utf-8") as f:
        f.write(sections)


def write_orientation_ditributions(
    run_folder, file="orientation_distributions.inp", ori_dist=None, default_ori=None
):
    s = ""
    for key, dist in ori_dist.items():
        s += f"*Orientation, name={key}\n{key}\n3, 0.\n"
        s += f"*Distribution, NAME={key}, location=ELEMENT, Table={key}_Table\n"
        # default value
        ori = default_ori[key.split("dist_")[-1]]["R"].T
        s += f", {ori[0][0]:0.6E}, {ori[1][0]:0.6E}, {ori[2][0]:0.6E}, {ori[0][1]:0.6E}, {ori[1][1]:0.6E}, {ori[2][1]:0.6E} \n"
        for id, ori in dist.items():
            s += f"{id:d}, {ori[0][0]:0.6E}, {ori[1][0]:0.6E}, {ori[2][0]:0.6E}, {ori[0][1]:0.6E}, {ori[1][1]:0.6E}, {ori[2][1]:0.6E} \n"

    # write string of input lines to file
    with open(run_folder / file, "w", encoding="utf-8") as f:
        f.write("".join(s))


def write_distribution_tables(
    run_folder, file="distribution_tables.inp", ori_dist=None
):
    """Tables need to be at the assembly level after the material definitions."""

    s = ""
    for key, dist in ori_dist.items():
        # print table
        s += f"*Distribution Table, NAME={key}_Table\ncoord3d, coord3d\n"

    # write string of input lines to file
    with open(run_folder / file, "w", encoding="utf-8") as f:
        f.write("".join(s))


def check_shell_element_normals(
    element_sets, elements, shell_sections, nodes, coordinate_systems, tol=90.0
):
    """For every element, check that the defined local coordinate system
    Z-direction is aligned with the shell normal within a set tolerance in degrees.
    Assumes quad or tria elements."""

    for name, section in shell_sections.items():
        elset = section["elset"]
        found_warnings = False
        for e in element_sets[elset]:
            # get xyz of nodes
            e_corner_xyz = np.array(
                [nodes[n]["global_xyz"] for n in elements[e]["nodes"]]
            )

            # determine approximate shell normal from nodes positions and order (normal defined by right-hand rule)
            shell_normal = np.cross(
                (e_corner_xyz[1, :] - e_corner_xyz[0, :]),
                (e_corner_xyz[2, :] - e_corner_xyz[1, :]),
            )

            # angle in degrees between the local orientation and the shell normal from the dot product
            ori_z = coordinate_systems[section["orientation"]]["R"][2]
            angle = np.rad2deg(
                np.arccos(
                    np.dot(shell_normal, ori_z)
                    / (np.linalg.norm(shell_normal) * np.linalg.norm(ori_z))
                )
            )

            # if angle exceeds tolearance value, issue warning
            if angle - tol > 0:
                warnings.warn(
                    f"Shell element {e} in set {elset} has normal that exceeds z-dir tolerance by {angle - tol} degrees!",
                    UserWarning,
                )
                found_warnings = True

        if not found_warnings:
            # check passed
            print(
                f"Shell element normals in set {elset} are within tol of local z-dir."
            )


def check_rts_path_angles(rts_plies, laminates, angle_tol=3):
    for lname, laminate in laminates.items():
        for ply in laminate["plies"]:
            if ply["name"] in rts_plies:
                # compare nominal UD material angle and calculated angle from steering path
                nominal_angle = ply["angle"]
                rts_path_angle = rts_plies[ply["name"]]["orientations"][
                    laminate["orientation"]
                ]["z_rot_angle"]

                # if angle exceeds tolearance value, issue warning
                if abs(rts_path_angle - nominal_angle) > angle_tol:
                    warnings.warn(
                        f"Control point path angle for base laminate {lname} and RTS ply {ply['name']} outside angle tolerance by {rts_path_angle - nominal_angle} degrees!",
                        UserWarning,
                    )
                else:
                    # check passed
                    print(
                        f"Control point path angle for base laminate {lname} and RTS ply {ply['name']} is within tol of {angle_tol} deg."
                    )


def plot_3D_scatter_properties(
    element_sets, elements, shell_sections, rts_discrete_reinforcements, view_projection
):
    if not view_projection:
        view_projection = {"elev": 90, "azim": 0, "roll": 90}

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    for name, section in shell_sections.items():
        p = np.empty([0, 3])
        for e in element_sets[section["elset"]]:
            p = np.append(p, [elements[e]["global_xyz"]], axis=0)
        ax.scatter(p[:, 0], p[:, 1], p[:, 2], label=name)

    if rts_discrete_reinforcements:
        for name, ply in rts_discrete_reinforcements.items():
            ctrl_points = np.array(
                [
                    np.array(ply["path_start"]) + p["d"] * ply["upath"]
                    for p in ply["control_pts"]
                ]
            )
            ax.plot(
                ctrl_points[:, 0],
                ctrl_points[:, 1],
                ctrl_points[:, 2],
                "-ko",
                label="control points",
            )

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_aspect("equal", adjustable="box")
    ax.view_init(**view_projection)
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")
    plt.show(block=True)


class Timer(object):
    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time.time()

    def __exit__(self, type, value, traceback):
        if self.name:
            print(
                "[%s]" % self.name,
            )
        print("Elapsed: %s" % (time.time() - self.tstart))


if __name__ == "__main__":

    # test inputs for arc_model.cae Model_1
    test_inputs = Path("../FE_model")
    run_folder = Path("./outputs")

    with open("./Study_1_1_flattened.json", "r") as f:
        rts_settings = json.load(f)

    PARAMETERS = [
        c["parameters"] for c in rts_settings["components"] if c["name"] == "abaqus"
    ][0]
    PART = "panel"

    for p in test_inputs.glob("**/*.inp"):
        copy2(p, run_folder / p.name)

    IGNORE_ELSETS = [
        # any non-design elsets, for example sets of sets that are only used for post-processing
    ]

    PLIES_EXCLUDED_FROM_THICKNESS_RATIO = [
        # e.g. outer woven pies that are not being designed and where thickness ratio constraints do not apply
    ]

    with Timer("setup"):
        setup(PART, PARAMETERS, run_folder)

    with Timer("compute"):
        compute(PART, PARAMETERS, run_folder)
