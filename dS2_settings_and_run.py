class ModelSettings:
    def __init__(self):
        # =====================================================================
        # General simulation values
        # =====================================================================
        # Start of the simulation period (yyyy, mm, dd, hh)
        self.sim_start      = 2010, 1, 1, 0
        # End of the simulation period, plus 1 dt (yyyy, mm, dd, hh)
        self.sim_end        = 2011, 1, 1, 0
        # Time step [h]
        self.dt             = 1

        # =====================================================================
        # File locations
        # =====================================================================
        # Location of the input files
        self.input_dir      = "input\\"
        # Location of the output folder
        self.output_dir     = "output\\"

        # =====================================================================
        # Catchment information (relative to input_dir)
        # =====================================================================
        # Location of (sub)basin by ID (matching info in routingInfo below), nan outside main basin [-]
        self.catchment      = "map_catchment.asc"
        # Distance of each pixel to main outlet [m]
        self.distance       = "map_distance_meters.asc"
        # Pickle object with 3 dictionaries: outletLoc, catchLoc, outletInfo
        self.routing_info   = "routingInfo.pkl"

        # =====================================================================
        # Forcing files (relative to input_dir)
        # =====================================================================
        # Numpy memmap files containing precipitation and evaporation data [mm h-1]
        self.precipitation  = "2010_tp.dat"
        self.evaporation    = "2010_etref.dat"
        # Start of the forcing data (yyyy, mm, dd, hh)
        self.data_start     = 2010, 1, 1, 0
        # End of the forcing data, plut 1 dt (yyyy, mm, dd, hh)
        self.data_end       = 2011, 1, 1, 0

        # =====================================================================
        # Initial conditions
        # =====================================================================
        # Initial snow storage and discharge value, can also be a NetCDF file
        self.init_Sstore    = 0
        self.init_Qsim      = 0.1

        # =====================================================================
        # Parameters
        # =====================================================================
        # Intersect of gQ function
        self.alpha          = -2.471
        # Slope of gQ function
        self.beta           = 0.861
        # Downward curvature of gQ function
        self.gamma          = -0.011
        # Correction factor on evaporation [-]
        self.eps            = 0.89
        # Average flow speed for routing [m s-1]
        self.tau            = 2

        # =====================================================================
        # Model options and settings
        # =====================================================================
        # Option to include snow processes in model
        self.SnowProcesses  = True
        # Numpy memmap file containing temperature [deg C]
        self.temperature    = "2010_t2m.dat"
        # Numpy memmap file containing radiation [W m-2]
        self.radiation      = "2010_ssrd.dat"
        # Critical temperature for snowfall and snowmelt [deg C]
        self.crit_temp      = 0
        # Degree day factor for snowmelt [mm day-1 deg C-1]
        self.ddf            = 2.00
        # Radiation day factor for snowmelt [mm day-1 (W/m2)-1]
        self.rdf            = 0.26

        # Option to include evaporation reduction
        self.EvapReduction  = True
        # Threshold below which no evaporation is allowed to occur
        self.Q_threshold    = 1e-4

        # =====================================================================
        # Reporting options
        # =====================================================================
        # Define which variables need to be stored as output, Qsim will always be stored as output
        self.output_Qsim    = True
        self.output_Sstore  = False # Only when SnowProcesses = True
        self.output_Smelt   = False # Only when SnowProcesses = True
        self.output_watbal  = False

        # =====================================================================
        # Numerics and stability settings
        # =====================================================================
        # Define approx maximum memory size that dS2 is allowed to use [MB]
        self.max_RAM        = 1024
        # Size of a single value in the memmap files: float32: 4, float64: 8 [bytes]
        self.size_one_value = 4
        # Q_t is not allowed to be smaller than Q_t-1 * LB
        self.LB             = 1e-3
        # Maximum allowed difference between gQ_t and gQ_t-1
        self.max_gQ_diff    = 2
        # Factor controlling how the time step is reduced
        self.dt_reduction   = 0.15
        # Minimum extra internal time steps
        self.min_extra_dt   = 5
        # Maximum extra internal time steps
        self.max_extra_dt   = 50

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from dS2_model import dS2
    # Read model settings
    Settings = ModelSettings()
    # Setup model
    model = dS2()
    model.model_setup_class(Settings)

    # Change model parameters using the functions below, or modify them in the settings.cfg file
    model.change_param("alpha", -3.471)
    model.change_param("beta", 0.561)
    model.change_param("gamma", -0.11)
    model.change_param("eps", 1.5)
    model.change_param("tau", 1.5)

    # Run model and rout water
    model.generate_runoff(progress=True)
    model.routing()

    Qsim = model.Qrout["1.0"][:len(model.sim_period)]

    # Plot results
    fig = plt.figure("Compare discharge", clear=True, tight_layout=True)
    ax = fig.subplots(1)
    ax.plot(Qsim, label="Simulated")
    ax.legend()