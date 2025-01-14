""" 
Execute this script from the cmd line with:
abaqus python extract_constraints_data_from_odb.py
"""
import sys
import json
from datetime import datetime

from odbAccess import *
from abaqusConstants import *
from odbMaterial import *
from odbSection import *
import numpy as np
from abaqusConstants import *

# read data from the settings.json file
with open("dapta_inputs.json", "r") as f:
    dapta_inputs = json.load(f)
open_mdao_in = [
        d for d in dapta_inputs["components"] if d["name"] == "open-mdao"
    ][0]
abaqus_in = [d for d in dapta_inputs["components"] if d["name"] == "abaqus"][0]
odb_path = str(abaqus_in["parameters"]["analysis_file"].replace(".inp",".odb"))

print odb_path

# set outputs hierarchy for lookup
OUTPUT_KEY_CHAIN = {
    "MASS":("step","instance","elset"),
    "U":("step","instance","node","dof"),
    "T-RATIO-MIN":("instance","elset","plyNo"),
    "T-RATIO-MAX":("instance","elset","plyNo"),
    "TSAIW":("step","instance","elset","plyNo"),
    "E":("step","instance","elset","plyNo","straintype")
    }

# sets initial values of the OUTPUTS 
OUTPUT_INIT = {"MASS":0.0, "U": 0.0, "TSAIW":0.0, "E":0.0}

# initialise the output object that will be written to file for further processing by the CLI
EXTRACTED_OUTPUTS = {}

# parse the abaqus outputs object into a local dictionary
# EXAMPLE parsed object
# {u'MASS': {u'Step-1': {u'PART-1-1': {u'ALL': 0.0, u'STEERED': 0.0}}}}
def create_dict(value, kchain, ref_dict):
    if len(kchain) == 0:
        return value
    elif len(kchain) == 1:
        if kchain[0] in ref_dict:
           ref_dict[kchain[0]].update(create_dict(value, {}, ref_dict[kchain[0]]))
           return ref_dict
        return {kchain[0]: create_dict(value, {}, {})}
    elif len(kchain)>=2:
        if kchain[0] in ref_dict:
           ref_dict[kchain[0]].update(create_dict(value, kchain[1:], ref_dict[kchain[0]]))
           return ref_dict
        return {kchain[0]: create_dict(value, kchain[1:], {})}
    
for key in abaqus_in["outputs"].keys():    
    kchain = key.split(".")

    output_n= kchain[0]
    
    if not output_n in OUTPUT_KEY_CHAIN:
        raise ValueError("Output type is not defined: " + str(output_n))  
    elif output_n == 'T-RATIO-MIN' or output_n == 'T-RATIO-MAX':
        pass 
    else:
        EXTRACTED_OUTPUTS.update(create_dict(OUTPUT_INIT[output_n],kchain,EXTRACTED_OUTPUTS))

# print str(EXTRACTED_OUTPUTS).encode('ascii') # DEBUG PRINT

# load simulation output database
odb = openOdb(path=odb_path)

# # '-----------------------------------------------------------------------------'
def get_mass(odb):
    
    for stepName in odb.steps.keys():
        
        if not EXTRACTED_OUTPUTS[outputName].get(stepName):
            continue # go to next step

        for instance_name in EXTRACTED_OUTPUTS[outputName].get(stepName):
            for elset_name in EXTRACTED_OUTPUTS[outputName].get(stepName).get(instance_name):
                
                data_found_flag = 0
                elset_str = str("ElementSet " + instance_name + "." + elset_name)                    
                houtputs = odb.steps[stepName].historyRegions[elset_str].historyOutputs
                for fieldName in houtputs.keys():                    
                    if fieldName == outputName:
                        EXTRACTED_OUTPUTS[outputName][stepName][instance_name][elset_name] = float(houtputs[fieldName].data[0][1])
                        data_found_flag +=1

                if not data_found_flag:
                    raise ValueError("Cannot find output for key chain: " + outputName + stepName + instance_name + elset_name)
                    
def get_displacement(odb):
    
     for stepName in odb.steps.keys():
        
        if EXTRACTED_OUTPUTS[outputName].get(stepName):
            lastFrame = odb.steps[stepName].frames[-1]
        else:
            continue # go to next step
        
        for instance_name in EXTRACTED_OUTPUTS[outputName].get(stepName):
            for nid in [nid for nid in EXTRACTED_OUTPUTS[outputName].get(stepName).get(instance_name) if nid.startswith("N")]:
                
                # check if the node id exists 
                try:
                    int(nid[1:])
                except:
                    ValueError("Cannot parse node id (use format N<nid>, e.g. N10): " + nid)

                # print nid[1:], outputName # DEBUG PRINT
                disp = lastFrame.fieldOutputs[str(outputName)].getSubset(position=NODAL).values[int(nid[1:])-1].data

                for dof in EXTRACTED_OUTPUTS[outputName].get(stepName).get(instance_name).get(nid):
                    data_found_flag = 0
                    
                    if dof in ["1", "2", "3"]:
                        EXTRACTED_OUTPUTS[outputName][stepName][instance_name][nid][dof]=float(disp[int(dof)-1])
                        data_found_flag +=1
                    elif dof == "norm":
                        EXTRACTED_OUTPUTS[outputName][stepName][instance_name][nid][dof]=float(np.linalg.norm(disp))
                        data_found_flag +=1
                
                    if not data_found_flag:
                        raise ValueError("Cannot find output for key chain: " + outputName + stepName + instance_name + nid) 
            
def get_tsaiw(odb):
    
    for stepName in odb.steps.keys():
        
        if not EXTRACTED_OUTPUTS[outputName].get(stepName):
            continue # go to next step

        for instance_name in EXTRACTED_OUTPUTS[outputName].get(stepName):
            for elset_name in EXTRACTED_OUTPUTS[outputName].get(stepName).get(instance_name):
                for ply_number in EXTRACTED_OUTPUTS[outputName].get(stepName).get(instance_name).get(elset_name):
                                        
                    data_found_flag = 0

                    try:
                        ele_set = odb.rootAssembly.instances[str(instance_name)].elementSets[str(elset_name)]
                    except:
                        raise ValueError("Element set not found: "+elset_name)
                    
                    spts_per_ply = 3 # >=2
                    category = odb.rootAssembly.instances[str(instance_name)].elementSets[str(elset_name)].elements[0].sectionCategory
                                                            
                    lastFrame = odb.steps[stepName].frames[-1]

                    try:
                        f = lastFrame.fieldOutputs[str(outputName)]
                    except:
                        raise ValueError(outputName+" could not be found in the requested Abaqus field outputs.")

                    spt_indices = (int(ply_number)*(spts_per_ply), int(ply_number)*(spts_per_ply) + spts_per_ply-1) # bottom and top surfaces
                    spt_crit = 0.0

                    for spt_index in spt_indices:
                        try:
                            spt = category.sectionPoints[spt_index]  
                        except:
                            raise ValueError("Check the requested ply: "+elset_name+"."+ply_number)                  
                        fieldValues=f.getSubset(region=ele_set, position=INTEGRATION_POINT, sectionPoint=spt).values
                        # find the element with the highest value per ply 
                        values_list = [(v.data) for v in fieldValues]
                        try:
                            max_data = max(values_list)
                        except:
                            raise ValueError("Section point error - Check the defined points in the requested field outputs.")
                        ele_index =  [(v.elementLabel) for v in fieldValues if v.data == max_data][0]  

                        if max_data > spt_crit:
                            spt_crit = max_data
                            EXTRACTED_OUTPUTS[outputName][stepName][instance_name][elset_name][ply_number]={"max": float(max_data), "ele": int(ele_index), "sectionPoint": spt_index}
                            data_found_flag +=1
                            # print outputName

def get_strain(odb):
    
    for stepName in odb.steps.keys():
        
        if not EXTRACTED_OUTPUTS[outputName].get(stepName):
            continue # go to next step

        for instance_name in EXTRACTED_OUTPUTS[outputName].get(stepName):
            for elset_name in EXTRACTED_OUTPUTS[outputName].get(stepName).get(instance_name):
                for ply_number in EXTRACTED_OUTPUTS[outputName].get(stepName).get(instance_name).get(elset_name):
                                        
                    data_found_flag = 0

                    try:
                        ele_set = odb.rootAssembly.instances[str(instance_name)].elementSets[str(elset_name)]
                    except:
                        raise ValueError("Element set not found: "+elset_name)
                    
                    spts_per_ply = 3 # >=2
                    category = odb.rootAssembly.instances[str(instance_name)].elementSets[str(elset_name)].elements[0].sectionCategory
                                                            
                    lastFrame = odb.steps[stepName].frames[-1]

                    try:
                        f = lastFrame.fieldOutputs[str(outputName)]
                    except:
                        raise ValueError(outputName+" could not be found in the requested Abaqus field outputs.")

                    spt_indices = (int(ply_number)*(spts_per_ply), int(ply_number)*(spts_per_ply) + spts_per_ply-1) # bottom and top surfaces
                    spt_crit = 0.0

                    for spt_index in spt_indices:
                        try:
                            spt = category.sectionPoints[spt_index]  
                        except:
                            raise ValueError("Check the requested ply: "+elset_name+"."+ply_number)                  
                        
                        
                        for straintype in EXTRACTED_OUTPUTS[outputName].get(stepName).get(instance_name).get(elset_name).get(ply_number):
                            
                            if straintype in ["MAX_INPLANE_PRINCIPAL", "MIN_INPLANE_PRINCIPAL", "MAX_PRINCIPAL", "MID_PRINCIPAL", "MIN_PRINCIPAL"]:
                                fieldValues=f.getSubset(region=ele_set, position=INTEGRATION_POINT, sectionPoint=spt).getScalarField(invariant=SymbolicConstant(straintype),).values
                            elif straintype in ["E11", "E22", "E12"]:
                                fieldValues=f.getSubset(region=ele_set, position=INTEGRATION_POINT, sectionPoint=spt).getScalarField(componentLabel=str(straintype),).values
                            else:
                                raise ValueError("Strain type not recognized: "+straintype)
                            
                            # find the element with the highest value per ply 
                            values_list = [(v.data) for v in fieldValues]
                            try:
                                max_data = max(values_list)
                            except:
                                raise ValueError("Section point error - Check the defined points in the requested field outputs.")
                            ele_index =  [(v.elementLabel) for v in fieldValues if v.data == max_data][0]  

                            if max_data > spt_crit:
                                spt_crit = max_data
                                EXTRACTED_OUTPUTS[outputName][stepName][instance_name][elset_name][ply_number][straintype]={"max": float(max_data), "ele": int(ele_index), "sectionPoint": spt_index}
                                data_found_flag +=1
                                # print outputName

# use a function selector dictionary to avoid nested if statements
get_outputs = {
    "MASS": get_mass, 
    "U": get_displacement,
    "TSAIW": get_tsaiw,
    "E": get_strain
}

# read the outputs from the database into the local dictionary    
for outputName in EXTRACTED_OUTPUTS.keys():
    
    try: 
        get_outputs[outputName](odb) # 
    except Exception as e:
        raise ValueError(str(e))
            
print str(EXTRACTED_OUTPUTS) # DEBUG PRINT

# dump data to file
EXTRACTED_OUTPUTS["timestamp"] = datetime.now().strftime("%Y-%m-%d-%H_%M_%S")
with open("odb_output.json", "w") as f:
    json.dump(EXTRACTED_OUTPUTS, f)
    
print "GOT RESULTS FROM ODB."