from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import xarray as xr
import pickle, os

def readMap(fileName, flat = False):
    """ Read 2D map, stored as an ASCII file, without any headers. Set flat to
    True if a flat array is required."""
    data = np.loadtxt(fileName, skiprows = 0)

    if flat == True:
        data = data.flatten()
        data = data[~np.isnan(data)]
    return data

def read_memmap(self, fileName, shape = None):
    """ Read the memmap file, given a certain shape. If the shape is not
    provided, the function will try to estimate it based on the the latitude
    map."""
    # Estimate the shape if none is provided
    if shape == None:
        # Determine the number of timesteps present in file
        d_length = len(np.memmap(fileName, dtype = "float32")) / \
                        len(self.catchment[~np.isnan(self.catchment)])
        # Set the shape correctly
        shape = (int(d_length), len(self.catchment[~np.isnan(self.catchment)]))
    # Read the memmap file
    data = np.memmap(fileName, dtype = "float32", mode = "r", shape = shape)
    return data

def read_class(self, settings):
    """ Read and extract the model settings from the settings class """
    # =========================================================================
    # General simulation values
    # =========================================================================
    s_start             = settings.sim_start
    s_end               = settings.sim_end
    self.dt             = settings.dt
    self.sim_period     = pd.DatetimeIndex(
            np.arange(start = datetime(*s_start),
                      stop  = datetime(*s_end),
                      step  = timedelta(hours=self.dt)).astype(datetime))
    self.tsteps         = len(self.sim_period)

    # =========================================================================
    # File locations
    # =========================================================================
    self.indir          = settings.input_dir
    self.outdir         = settings.output_dir
    if not os.path.exists(self.outdir): os.mkdir(self.outdir)

    # =========================================================================
    # Catchment information
    # =========================================================================
    self.catchment      = readMap(self.indir + settings.catchment)
    self.shape          = (self.tsteps, len(self.catchment[~np.isnan(self.catchment)]))
    # Get distance information
    catch1D             = self.catchment.flatten()
    catch1D             = catch1D[~np.isnan(catch1D)]
    self.dist1D         = np.zeros(catch1D.shape)
    self.dist1D[self.dist1D == 0] = np.nan
    self.dist1D[catch1D >= 1] = readMap(self.indir + settings.distance, flat = True)
    # Read pickle object
    self.outletLoc, self.catchLoc, self.outletInfo = pickle.load(
            open(self.indir + settings.routing_info, "rb"))

    # =========================================================================
    # Forcing files
    # =========================================================================
    # Read precipitation and evaporation
    self.P_file         = self.indir + settings.precipitation
    self.ET_file        = self.indir + settings.evaporation
    # Timing data related to forcing files
    d_start             = settings.data_start
    d_end               = settings.data_end
    self.data_period    = pd.DatetimeIndex(
            np.arange(start = datetime(*d_start),
                      stop = datetime(*d_end),
                      step = timedelta(hours=self.dt)).astype(datetime))
    # =========================================================================
    # Initial conditions
    # =========================================================================
    for var in ["Qsim", "Sstore"]:
        # Try to read initial storage (first as value, than as map, than as NetCDF)
        try:
            dat         = float(getattr(settings, "init_{}".format(var)))
            setattr(self, "init_{}".format(var), dat)
        except:
            try:
                # print(self.indir + getattr(settings, var))
                dat = readMap(self.indir + getattr(settings, "init_{}".format(var)), flat=True)
                setattr(self, "init_{}".format(var), dat)
            except:
                fileLoc     = self.indir + getattr(settings, "init_{}".format(var))
                dat         = readNetCDF(fileLoc).loc[self.sim_period[0]]
                dat         = np.array(dat).flatten()
                setattr(self, "init_{}".format(var), dat[~np.isnan(dat)])

    # =========================================================================
    # Parameters
    # =========================================================================
    for var in ["alpha", "beta", "gamma", "eps"]:
        # Try to read as a single value or as a map
        dat             = getattr(settings, var)
        try:
            value       = float(dat)
        except:
            value       = readMap(self.indir + dat, flat = True)
        self.change_param(var, value)
    # Only a single value is accepted for the routing parameter
    self.change_param("tau", settings.tau)

    # =========================================================================
    # Model options and settings
    # =========================================================================
    # Snow processes
    self.SnowFLAG       = settings.SnowProcesses
    if self.SnowFLAG:
        # Read additional file locations
        self.T_file     = self.indir + settings.temperature
        self.Rs_file    = self.indir + settings.radiation
        # Set parameters
        self.change_param("T0", settings.crit_temp)
        self.change_param("ddf", settings.ddf)
        self.change_param("rdf", settings.rdf)

    # Evaporation reduction
    self.evap_reduction = settings.EvapReduction
    if self.evap_reduction:
        self.change_param("Qt", settings.Q_threshold)

    # =========================================================================
    # Reporting options
    # =========================================================================
    # Qsim is always true, as this is required for routing
    self.Qsim_flag      = True
    self.Qsim_flag_nc   = settings.output_Qsim
    self.Sstore_flag    = settings.output_Sstore
    self.Smelt_flag     = settings.output_Smelt
    self.WatBal_flag    = settings.output_watbal

    # Keep track of how many files are simultaneously in memory
    # P (forcing), E (evaporation), and Q (output) are always in memory
    file_count          = 3

    # Check how many output files are simultaneously in memory
    for var in ["Qsim_flag", "Smelt_flag", "Sstore_flag", "WatBal_flag"]:
         if getattr(self, var):
             file_count += 1

    # =========================================================================
    # Numerics and stability settings
    # =========================================================================
    # Check the number of forcing files, dependent on the model options
    # Snow uses the 2 forcing data (T and Rg), so thats why an elif
    if self.SnowFLAG:
        file_count     += 2
    # Convert max RAM from megabytes to bytes
    max_RAM             = settings.max_RAM * 1e6
    size_one_value      = settings.size_one_value
    # Determine data type based on the size of one value
    if size_one_value == 4:
        self.dtype      = np.float32
    elif size_one_value == 8:
        self.dtype      = np.float64
    else:
        raise ValueError("Value of {:.1f} for size_one_value is not valid".format(size_one_value))
    # Determine the chunk size
    max_size            = max_RAM//file_count
    self.chunk_size     = int(max_size / (size_one_value * self.shape[1]))
    # Read other parameters
    self.LB             = settings.LB
    self.max_gQ_difference = settings.max_gQ_diff
    self.dt_reduction   = settings.dt_reduction
    self.min_extra_dt   = settings.min_extra_dt
    self.max_extra_dt   = settings.max_extra_dt


def readNetCDF(fileName):
    """ Read netCDF files form disk and return the variable that is inside """
    data = xr.open_dataset(fileName, autoclose = True)
    var = list(data.data_vars.keys())[0]
    data = data.rename({var: 'val'})
    return data.val