""" Optimisation Component classes for OpenMDAO and associated utilities."""

import numpy as np
import openmdao.api as om  # type: ignore

from component_api2 import call_compute, call_setup
import traceback


class OMexplicitComp(om.ExplicitComponent):
    """standard component that follows the OM conventions"""

    def __init__(self, compname, fd_step, has_compute_partials=True):
        super().__init__()
        self.compname = compname
        self.get_grads = True
        self.iter = 0
        self.partial_dict = None
        self.fd_step = fd_step
        self._has_compute_partials = has_compute_partials  # overrides parent class attribute with user defined parameter

    def setup(self):
        message = {"component": self.compname}
        _, component_dict = call_setup(message)
        inputs = component_dict["input_data"]["design"]
        outputs = component_dict["output_data"]["design"]

        # initialise the inputs
        if inputs:
            for variable in inputs:
                self.add_input(variable.replace(".", ":"), val=inputs[variable])

        # initialise the outputs
        if outputs:
            for variable in outputs:
                self.add_output(variable.replace(".", ":"), val=outputs[variable])

    def setup_partials(self):
        # Get the component partial derivative information

        message = {"component": self.compname}
        _, component_dict = call_setup(message)

        if "partials" in component_dict and component_dict["partials"]:
            self.partial_dict = component_dict["partials"]
            for resp, vars in self.partial_dict.items():
                for var, vals in vars.items():
                    self.declare_partials(
                        resp.replace(".", ":"), var.replace(".", ":"), **vals
                    )
        else:
            # calculate all paritials using finite differencing
            self.declare_partials("*", "*", method="fd", step=self.fd_step)

    def compute(self, inputs, outputs, discrete_inputs=None, discrete_outputs=None):
        print("Calling compute.")
        # calculate the outputs
        # Note: transform all np.ndarrays into nested lists to allow formatting to json
        input_dict = {"design": reformat_inputs(inputs._copy_views())}

        message = {
            "component": self.compname,
            "inputs": input_dict,
            "get_grads": True,
            "get_outputs": True,
        }
        print("message: \n", str(message))

        try:
            _, data = call_compute(message)
            if not "outputs" in data:
                raise ValueError(f"Error: Compute output missing - output was: {data}.")
        except Exception as e:
            print(f"Compute of {self.compname} failed, input data was: {str(message)}")
            tb = traceback.format_exc()
            print(tb)
            raise ValueError(
                f"OM Explicit component {self.compname} compute error: " + tb
            )

        val_outputs = data["outputs"]["design"]

        # OpenMDAO doesn't like the outputs dictionary to be overwritten, so
        # assign individual outputs one at a time instead
        for output in outputs:
            if "." in output:
                name = output.split(".")[1]  # remove the OpenMDAO component prefix
            else:
                name = output
            name = name.replace(":", ".")
            outputs[output] = val_outputs[name]

    def compute_partials(self, inputs, J):
        """Jacobian of partial derivatives."""

        print("Calling compute_partials.")

        self.iter += 1
        #input_dict = reformat_inputs(inputs._copy_views())
        input_dict = {"design": reformat_inputs(inputs._copy_views())}

        message = {
            "component": self.compname,
            "inputs": input_dict,
            "get_grads": True,
            "get_outputs": True,
        }
        print("message: \n", str(message))

        try:
            _, data = call_compute(message)
            if not "partials" in data:
                raise ValueError(
                    f"Error: Compute partial derivatives missing - output was: {data}."
                )
            self.partial_dict = data["partials"]
        except Exception as e:
            print(
                f"Compute partials of {self.compname} failed, input data was: {str(message)}"
            )
            tb = traceback.format_exc()
            print(tb)
            raise ValueError(
                f"OM Explicit component {self.compname} compute error: " + tb
            )

        if self.partial_dict:
            for resp, vars in self.partial_dict.items():
                for var, vals in vars.items():
                    if "val" in vals:
                        J[resp.replace(".", ":"), var.replace(".", ":")] = vals["val"]
            # print(dict(J))
        else:
            raise ValueError(f"Component {self.compname} has no Jacobian defined.")


class OMimplicitComp(om.ImplicitComponent):
    """standard implicit component that follows the OM conventions"""

    def __init__(self, compname, fd_step, has_compute_partials=False):
        super().__init__()
        self.compname = compname
        self.get_grads = True
        self.iter = 0
        self.partial_dict = None
        self.fd_step = fd_step

    def setup(self):
        message = {"component": self.compname}
        _, component_dict = call_setup(message)
        inputs = component_dict["input_data"]["design"]
        outputs = component_dict["output_data"]["design"]

        # initialise the inputs
        if inputs:
            for variable in inputs:
                self.add_input(variable.replace(".", ":"), val=inputs[variable])

        # initialise the outputs
        if outputs:
            for variable in outputs:
                self.add_output(variable.replace(".", ":"), val=outputs[variable])

    def setup_partials(self):
        message = {"component": self.compname}
        _, component_dict = call_setup(message)

        if "partials" in component_dict and component_dict["partials"]:
            self.partial_dict = component_dict["partials"]
            for resp, vars in self.partial_dict.items():
                for var, vals in vars.items():
                    self.declare_partials(
                        resp.replace(".", ":"), var.replace(".", ":"), **vals
                    )
        else:
            # calculate all paritials using finite differencing
            self.declare_partials("*", "*", method="fd", step=self.fd_step)

    def apply_nonlinear(self, inputs, outputs, residuals):
        input_dict = {"design": reformat_inputs(inputs._copy_views())}

        message = {
            "component": self.compname,
            "inputs": input_dict,
            "get_grads": True,
            "get_outputs": True,
        }
        print("message: \n", str(message))

        try:
            _, data = call_compute(message)
            if not "outputs" in data:
                raise ValueError(f"Error: Compute output missing - output was: {data}.")
        except Exception as e:
            print(f"Compute of {self.compname} failed, input data was: {str(message)}")
            tb = traceback.format_exc()
            print(tb)
            raise ValueError(
                f"OM Explicit component {self.compname} compute error: " + tb
            )

        val_outputs = data["outputs"]["design"]

        for output in outputs:
            residuals[output] = val_outputs[output.replace(":", ".")]


def reformat_inputs(inputs):
    input_dict = inputs
    for key in [*input_dict.keys()]:
        new_key = key.split(".")[-1].replace(":", ".")
        input_dict[new_key] = input_dict.pop(key)
        if isinstance(input_dict[new_key], np.ndarray):
            input_dict[new_key] = input_dict[new_key].tolist()
    return input_dict
