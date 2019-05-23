import numpy as np
from funcSolver import flexible_RK4

def snowmelt_potential(T, Rg, T0, ddf, fact):
    """
    Calculate the potential snowmelt for each chunk of forcing data, where the
    melt is not corrected for the current snow storage (this is done in the
    dynamic loop of the model). Snowmelt is the cumulative snowmelt during each
    time step, and not an intensitiy. Function is based on [1]_.

    Parameters
    ----------
    T : numpy array
        Array containing the air temperature of each model cell for that specific
        timestep
    Rg : numpy array
        Array containing the global radiation of each model cell for that
        specific timestep
    T0 : float
        Critical temperature, above which snow will melt
    ddf : float
        Degree day factor
    fact : float
        Conversion factor from W/m2 to mm/hour (0.026 as default)

    References
    ----------
    .. [1]  Kustas, W. P., Rango, A., and Uijlenhoet, R.: A Simple Energy Budget
    Algorithm for the Snowmelt Runoff Model,Water Resources Research, 30,
    1515–1527, https://doi.org/10.1029/94WR00152, 1994.
    """
    melt = ddf * np.fmax(0, T - T0) + fact * Rg
    return melt

def evaporation_reduction(QwithET, Q_init, P, alpha, beta, gamma, Q_threshold, LB, dt):
    ''' Corrects the pixels where Q < Q_threshold to have no evaporation. This
    way, negative discharge values are prevented. The function returns the new
    discharge timeseries, and the indices where evaporation reduction has
    occurred. If no evaporation reduction has occurred, the indiced variable
    will return None.'''
    # Extract the lowest value from the discharge estimate for the next time step
    Qmin = np.min(QwithET)
    # variable to store the indices where reduction has taken place
    idx = None
    # only perform evaporation reduction if the discharge threshold is crossed
    if Qmin < Q_threshold:
        # find the indices where the reduction needs to occur
        idx = np.where(QwithET <= Q_threshold)
        # extract the correct initial discharge and precipitation values
        Q_slice = Q_init[idx]
        precip_slice = P[idx]
        alpha = alpha[idx]
        beta = beta[idx]
        gamma = gamma[idx]
        boundary = LB
        # calculate the new discharge estimate without evaporation
        Qnew, extra_dt = flexible_RK4(P=precip_slice, ET=0, Q=Q_slice, alpha=alpha, beta=beta, gamma=gamma, LB=boundary, dt=dt)
        # update the values in the vector with the new values
        QwithET[idx] = Qnew
    return QwithET, idx