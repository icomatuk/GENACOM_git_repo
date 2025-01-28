""" Local component API. """

#import abaqus
#import link_variables
import paraboloid



def call_setup(msg):
    # if msg["component"] == "abaqus":
    #     return ("api completed call_setup", abaqus.setup(msg))
    # elif msg["component"] == "link_variables":
    #     return ("api completed call_setup", link_variables.setup(msg))
    if msg["component"] == "paraboloid":
        return ("api completed call_setup", paraboloid.setup(msg))    
    else:
        raise ValueError("Component not implemented.")


def call_compute(msg):
    # if msg["component"] == "abaqus":
    #     return ("api completed call_compute", abaqus.compute(msg))
    # elif msg["component"] == "link_variables":
    #     return ("api completed call_compute", link_variables.compute(msg))
    if msg["component"] == "paraboloid":
        return ("api completed call_compute", paraboloid.compute(msg))    
    else:
        raise ValueError("Component not implemented.")
