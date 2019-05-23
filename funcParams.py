import numpy as np
import pandas as pd
from datetime import datetime, timedelta

def set_parameter(self, name, value):
    '''Allows user to changes parameter value, and where the function will take
    care to correct the parameter value to the correct units. '''
    # List of all allowed parameters
    all_parameters = (
        "alpha", "beta", "gamma", "eps", "tau",
        "T0", "ddf", "rdf", "Qt",
        "sim_start", "sim_end", "dt", "Sstore_init", "Qsim_init",
        "max_RAM", "size_one_value", "LB",  "max_gQ_difference",
        "dt_reduction", "min_extra_dt", "max_extra_dt",
        "Qsim_flag", "Smelt_flag", "Sstore_flag", "watbal_flag")

    # Some parameters are not allowed to be changed on the fly, since they
    # require new input data
    not_allowed = ("dt", "max_RAM", "size_one_value")

    # Check whether it is allowed to change this parameter
    if name not in all_parameters:
        raise NameError(("Variable {} is not recognised as a parameter, "
                        "please try again with a different name. "
                        "Select one from this list: \n{}").format(name,
                                                              all_parameters))
    elif name in not_allowed:
        raise Exception(("Variable {} is not allowed to be changed on the fly "
                        "since it requires new input data or reloading of "
                        "current data in memory").format(name))
    # Update the start of the simulation
    if name == "sim_start":
        self.sim_period = pd.DatetimeIndex(
            np.arange(start = datetime(*value),
                      stop = self.sim_period[-1].to_pydatetime()
                                 + timedelta(hours=self.dt),
                      step = timedelta(hours=self.dt)).astype(datetime))
        self.tsteps = len(self.sim_period)
        self.shape[0] = self.tsteps
        # self.shape = (self.tsteps, len(self.Lat))
    # Update the end of the simulation
    elif name == "sim_end":
        self.sim_period = pd.DatetimeIndex(
            np.arange(start = self.sim_period[0].to_pydatetime(),
                      stop = datetime(*value),
                      step = timedelta(hours=self.dt)).astype(datetime))
        self.tsteps = len(self.sim_period)
        self.shape[0] = self.tsteps
        # self.shape = (self.tsteps, len(self.Lat))
    # Update the intinal conditions
    elif name == "Sstore_init":
        self.init_Sstore = value
    elif name == "Qsim_init":
        self.init_Qsim = value
    # Update t_lag parameter
    elif name == "tau":
        # Assume value is given in m/s
        lag_dt_m = 1/(value * (3600 * self.dt))
        self.tau = lag_dt_m
    elif name == "Qsim_flag":
        self.Qsim_flag_nc = value
    elif name == "ddf":
        # Assume value is given as a degree-day-factor
        self.ddf = value * self.dt/24
    elif name == "rdf":
        # Asuume value is given as a radiation-day-factor
        self.rdf = value * self.dt/24
    # If no different operations are required, just change the value
    else:
        setattr(self, name, value)
