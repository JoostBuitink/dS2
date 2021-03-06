import numpy as np
import pandas as pd

def with_diffusion(self, Qsim, main_ID=1.0, trackwater=False):

    # Extract all outlet values from dictionary
    allOutlets = list(map(float, self.outletLoc.keys()))
    # Create empty dictionary to store data in
    datalength = sum([int(item[1]*self.tau) for item in self.outletInfo.values()])
    Qtotal = {str(x): np.repeat(0., len(Qsim) + datalength) for x in allOutlets}

    # Loop through all outlets
    for outlet in allOutlets:
        # Convert to string
        outlet = str(outlet)

        # Find the index of the outlet
        outlet_ind = self.outletLoc[outlet]
        # Find the distance corresponding to this outlet
        outlet_dist = self.dist1D[outlet_ind]

        # Loop through all indices related to that outlet
        for ind in self.catchLoc[outlet]:
            ind = ind[0]
            # Extract simulated values and the corresponding timelag
            orig_values = Qsim[:,ind]
            # tmplag = max(int(self.dist1D[ind] * self.tau) - baseline_dist, 0)
            pixel_to_outlet = self.dist1D[ind] - outlet_dist
            #TODO: Fix errors where the timelag is less than the pixels (wrong input data)
            time_lag = max(0, int(pixel_to_outlet * self.tau))

            # Apply moving average
            window = int(time_lag * self.lag_to_window)
            values = apply_moving_window(window=window, values=orig_values)

            # Shift the simulated discharge with a timelag
            Qtotal[outlet][time_lag:time_lag+len(values)] += values

            # To other downstream outlets
            if outlet != str(main_ID):

                dwnID = self.outletInfo[outlet][0]
                # Continue adding discharge to downstream outlets until main outlet is reached
                while dwnID != "nan":
                    # Find the index of the downstream outlet
                    down_ind = self.outletLoc[dwnID]
                    # Find the distance corresponding to this outlet
                    down_dist = self.dist1D[down_ind]

                    # Calculate distance from pixel to outlet, and determine time lag
                    tmp_pixel_to_outlet = self.dist1D[ind] - down_dist
                    pixel_to_outlet = tmp_pixel_to_outlet
                    delay = int(pixel_to_outlet * self.tau)

                    # Apply moving average
                    window = int(delay * self.lag_to_window)
                    values = apply_moving_window(window=window, values = orig_values)

                    # Add values to the time series of the downstream outlet
                    Qtotal[dwnID][delay:delay+len(values)] += values

                    # Tracking of water
                    if dwnID == str(main_ID) and trackwater:
                        self.Qorigin_tmp.loc[delay:delay+len(values)-1, outlet] += values

                    # Update to next downstream outlet
                    dwnID = self.outletInfo[dwnID][0]
            else:
                if trackwater:
                    self.Qorigin_tmp.loc[time_lag:time_lag+len(values)-1, outlet] += values

    return Qtotal

def moving_average(values, window):
    ret = np.cumsum(values)
    ret[window:] = ret[window:] - ret[:-window]
    return ret[window - 1:] / window

def apply_moving_window(values, window):
    # Set zeroes, so only places where the window covers enough data gets values
    res = np.zeros(values.shape)
    if window > 3:
        ret = moving_average(values=values, window=window)
        # Find the correct indices to place the values in (the middle)
        start = window//2
        stop = -start + 1 if window%2 == 0 else -start
        # Copy original values, and replace with the averaged values
        res[start:stop] = ret
        res[values == 0] = 0
    # No need to define end, only causes errors
    elif window == 2:
        ret = moving_average(values=values, window=window)
        # Copy original values, and replace with the averaged values
        res[1:] = ret
        res[values == 0] = 0
    # If window is only 1, just copy the values
    else:
        res = values

    return res

def subcatch_to_outlet(self, Qsim, dist1D, outletLoc, catchLoc, t_lag):
    """
    Returns a dictionary with the runoff generated in each subbasin, routed to
    the outlet of the corresponding outlet. The outletIDs are used as keys.

    Parameters
    ----------
    Qsim : numpy ndarray
        A 2 dimension array containing the discharge generated per grid cell,
        per timestep, where each timestep is represented as a single row
    dist1D : numpy array
        An array containing the distance per cell towards the main outlet of
        the catchment (measured in meters)
    outletLoc : dict
        Dictionary containing the outlet IDs as keys and the indices related to
        these outlets as values. The indices are based on the shape and dimensions
        of the catchment map
    catchLoc : dict
        Dictionary containing the outlet IDs as keys and a list of indices
        where the subcatchment ID == outlet ID
    t_lag : float
        Value describing the time it takes for water to move 1 meter downstream
        in number of timesteps

    Returns
    -------
    Qrout : dict
        Dictionary containg the routed discharge of each subbasin to the
        corresponding outlet

    See Also
    --------
    :func:`through_allOutlets`: rout water through all subbasins

    Examples
    --------
    >>> Qsim = np.array([[7e0, 0,   0,   0, 0,   1e2, 2e1],
    ...                  [0,   0,   1e2, 0, 1e2, 0,   0],
    ...                  [0,   2e2, 0,   0, 0,   0,   0]])
    >>> dist1D = np.array([0, 1000, 2000, 3000, 2000, 3000, 4000])
    >>> outletLoc = {
    ...     '1.0': np.array([0]),
    ...     '2.0': np.array([6]),
    ...     '3.0': np.array([1])}
    >>> catchLoc = {
    ...     '1.0': np.array([[0]]),
    ...     '2.0': np.array([[6]]),
    ...     '3.0': np.array([[1], [2], [3], [4], [5]])}
    >>> t_lag = 0.001
    >>> result = subcatch_to_outlet(Qsim, dist1D, outletLoc, catchLoc, t_lag)
    >>> Qrout = result
    >>> Qrout
    {'1.0': array([7., 0., 0., 0., 0., 0., 0.]),
     '2.0': array([20., 0., 0., 0., 0., 0., 0.]),
     '3.0': array([0., 0., 500., 0., 0., 0., 0.])}
    """

    # Extract all outlet values from dictionary
    allOutlets = list(map(float, outletLoc.keys()))
    # Create empty dictionary to store data in
    Qrout = {}
    # Loop through all outlets
    for outlet in allOutlets:
        # Create empty vector to store data in
        Qrout[str(outlet)] = np.repeat(0., len(Qsim) + int(np.nanmax(dist1D) * t_lag))

        # Find the index of the outlet
        outlet_ind = outletLoc[str(outlet)]
        # Find the timelage corresponding to this outlet
        baseline_dist = int(dist1D[outlet_ind] * t_lag)
        # Loop through all indices related to that outlet
        for ind in catchLoc[str(outlet)]:
            ind = ind[0]
            # Extract simulated values and the corresponding timelag
            values = Qsim[:,ind]
            tmplag = max(int(dist1D[ind] * t_lag) - baseline_dist, 0)

            # Shift the simulated discharge with a timelag
            Qrout[str(outlet)][int(tmplag):int(tmplag)+len(values)] += values

    return Qrout

def through_allOutlets(self, outletINFO, Qrouted, catchLoc, t_lag, main_ID=1.0, trackwater=False):
    """
    Return the routed discharge per subbasin, taking a timelag into account.
    All values should be in mm/timestep

    Parameters
    ----------
    outletINFO : dict
        Dictionary containing the outlet indentifiers as keys, and a list of
        two items per outlet with (1) the ID of the downstream outlet and (2)
        the distance towards this outlet
    Qrouted : dict
        Dictionary containing the routed discharge per subbasin; where only water
        in the relevant subcatchment has been routed to this outlet
    catchLoc : dict
        Dictionary containing the outlet IDs as keys and a list of indices
        where the subcatchment ID == outlet ID
    t_lag : float
        Value describing the time it takes for water to move 1 meter downstream
        in number of timesteps

    Returns
    -------
    Qtotal : dict
        Dictionary containing the fully routed discharge, with both timelag and
        nested catchments taken into account in mm/timestep.

    See Also
    --------
    :func:`subcatch_to_outlet`: retrieve the Qrout and allIndices values

    Examples
    --------
    Chosen unrealistic large values Qrouted values to keep track of where the
    water is orginated from

    >>> outletINFO = {"1.0": ["nan", 0],
    ...               "2.0": ["3.0", 3000],
    ...               "3.0": ["1.0", 1000]}
    >>> Qrouted = {"1.0": [7e0, 0., 0.],
    ...            "2.0": [2e1, 0., 0.],
    ...            "3.0": [0., 0., 5e2]}
    >>> catchLoc = {
    ...     '1.0': np.array([[0]]),
    ...     '2.0': np.array([[6]]),
    ...     '3.0': np.array([[1], [2], [3], [4], [5]])}
    >>> t_lag = 0.001
    >>> Qtotal = through_allOutlets(outletINFO, Qrouted, catchLoc, t_lag)
    >>> Qtotal
    {'1.0': array([1., 0., 0., 71.42857143, 2.85714286, 0., 0.]),
     '2.0': array([20., 0., 0., 0., 0., 0., 0.]),
     '3.0': array([0., 0., 83.333, 3.333, 0., 0., 0.])}

    """
    # Dictionary to keep track of how many times water was added to each outlet
    count = {x: len(catchLoc[x]) for x in Qrouted}

    # Dictionary to store routed water in
    datalength = sum([int(item[1]*t_lag) for item in outletINFO.values()])
    Qtotal = {x: np.insert(np.repeat(0.,datalength), 0, Qrouted[x]) for x in Qrouted}

    # Loop through all outlets
    for outletID in outletINFO:
        # Exclude the main outlet
        if outletID != str(main_ID):
            # Get the ID of the downstream catchment (and thus outlet)
            dwnID = outletINFO[outletID][0]
            values = Qrouted[outletID]

            # Get the timelag between the current outlet and the downstream outlet
            delay = int(outletINFO[outletID][1] * t_lag)

            # Only if there is a downstream outlet
            while dwnID != "nan":
                # Add the flow of the outletID to the downstream outlet
                Qtotal[dwnID][delay:delay+len(values)] += values

                if dwnID == str(main_ID) and trackwater:
                    self.Qorigin_tmp[outletID][delay:delay+len(values)] = values

                # Keep track of how many times this addition has taken place
                count[dwnID] += len(catchLoc[outletID])
                # Get the delay to the downstream outlet and set the downstreamID
                # to the next downstream outlet
                delay += int(outletINFO[dwnID][1] * t_lag)
                dwnID = outletINFO[dwnID][0]
        else:
            if trackwater:
                values = Qrouted[outletID]
                delay = 0
                self.Qorigin_tmp[outletID][delay:delay+len(values)] = values

    Qtotal = {x: Qtotal[x]/count[x] for x in Qtotal}

    return Qtotal

def multiOutlet_routing(self, Qsim, main_ID=1.0, trackwater=False):
    """
    Returns the routed discharge for each individual outlet. This function is
    a wrapper for three other functions: subcatch_to_outlet, get_outletINFO,
    through_allOutlets.

    Parameters
    ----------
    Qsim : numpy ndarray
        A 2D array with a map of generated discharge per gridcell, with a row
        of values per timestep

    Returns
    -------
    Qmulti : dict
        Dictionary with discharge timeseries per outlet, where the ID of the
        outlet is used as key.

    See Also
    --------
    :func:`subcatch_to_outlet`: retrieve the Qrout and allIndices values
    :func:`through_allOutlets`: rout water through all subbasins

    Examples
    --------
    See the individual functions for examples
    """

    Qrouted = subcatch_to_outlet(self, Qsim=Qsim,
                                 dist1D=self.dist1D,
                                 outletLoc=self.outletLoc,
                                 catchLoc=self.catchLoc,
                                 t_lag=self.tau)

    Qmulti = through_allOutlets(self, outletINFO=self.outletInfo,
                                Qrouted=Qrouted,
                                catchLoc=self.catchLoc,
                                t_lag=self.tau,
                                main_ID=main_ID,
                                trackwater=trackwater)

    return Qmulti