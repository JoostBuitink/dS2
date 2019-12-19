import numpy as np

def readMap(fileName, flat = False):
    """ Read 2D map, stored as an ASCII file, without any headers. Set flat to
    True if a flat array is required."""
    try:
        data = np.loadtxt(fileName, skiprows = 0)
    except ValueError:
        data = np.loadtxt(fileName, skiprows = 6)

    if flat == True:
        data = data.flatten()
        data = data[~np.isnan(data)]
    return data

def correct_NA(data, NA_value):
    data[data == NA_value] = np.nan
    return data

def flat(data):
    data = data.flatten()
    data = data[~np.isnan(data)]
    return data

def get_neighbours(p, exclude_p=True, shape=None):
    """
    Returns the indices of cells neighbouring to cell with index P in a 2d array

    Parameters
    ----------
    p : numpy array
        Array with the index of the cell of interest
    exclude p : boolean flag
        Whether p should be ignored in the output
    shape : tuple
        Tuple with the shape of the array where the neighbours required

    Returns
    -------
    neighbours : array
        Array with the indices of the cell s neighbouring to p

    Examples
    --------
    >>> pp = np.r_[4, 5]
    >>> shape = (6, 6)
    >>> neighbours = get_neighbours(p, shape=shape)
    >>> x = np.zeros(shape, int)
    >>> x[tuple(neighbours.T)] = 1
    >>> x[tuple(p)] = 2
    >>> x
    [[0 0 0 0 0 0],
     [0 0 0 0 0 0],
     [0 0 0 0 0 0],
     [0 0 0 0 1 1],
     [0 0 0 0 1 2],
     [0 0 0 0 1 1]]

    """
    ndim = len(p)
    # generate an (m, ndims) array containing all combinations of 0, 1, 2
    offset_idx = np.indices((3,) * ndim).reshape(ndim, -1).T
    # use these to index into np.array([-1, 0, 1]) to get offsets
    offsets = np.r_[-1, 0, 1].take(offset_idx)
    # optional: exclude offsets of 0, 0, ..., 0 (i.e. p itself)
    if exclude_p:
        offsets = offsets[np.any(offsets, 1)]
    # Apply offsets to p
    neighbours = p + offsets
    # optional: exclude out-of-bounds indices
    if shape is not None:
        valid = np.all((neighbours < np.array(shape)) & (neighbours >= 0), axis=1)
        neighbours = neighbours[valid]
    return neighbours

def subcatch_to_outlet(distMap, catchmentMap, outletMap, elevMap):
    """
    Returns (1) a dictionary per outlet, the water in the corresponding subbasing
    routed to the outlet, (2) a dictionary with the indices of the outlets
    according to the shape of catchmentMap and (3) a dictionary with the indices
    of the subcatchments corresponding with the outlet ID according to the
    shape of catchmentMap.

    Parameters
    ----------
    distMap : numpy ndarray
        A 2 dimensional array containing the distance per cell towards the main
        outlet of the catchment (counted in number of timesteps)
    catchmentMap : numpy ndarray
        A 2 dimensional array containing the catchment identifiers, corresponding
        to the outlet identifiers
    outletMap : numpy ndarray
        A 2 dimensional array containing the location of the outlet, indicated
        with their identifiers. All cells without outlet are indicated with a zero

    Returns
    -------
    outletIndices : dict
        Dictionary containing the outlet IDs as keys and the indices related to
        these outlets as values. The indices are based on the shape and dimensions
        of the catchment map
    allIndices : dict
        Dictionary containing the outlet IDs as keys and a list of indices
        where the subcatchment ID == outlet ID

    See Also
    --------
    :func:`get_outletINFO`: retrieve the information about the outlets
    :func:`through_allOutlets`: rout water through all subbasins

    Examples
    --------
    >>> distMap = np.array([[np.nan, 0, np.nan],
    ...                     [1,      2, 3],
    ...                     [2,      3, 4]])
    >>> catchmentMap = np.array([[np.nan, 1., np.nan],
    ...                          [3.,     3., 3.],
    ...                          [3.,     3. ,2.]])
    >>> outletMap = np.array([[0, 1, 0],
    ...                       [3, 0, 0],
    ...                       [0, 0, 2]])
    >>> result = subcatch_to_outlet(Qsim, distMap,
    ...                        catchmentMap, outletMap)
    >>> outletIndices, allIndices, outletINFO = result
    >>> outletIndices
    {'1.0': array([0, 1]),
     '2.0': array([2, 2]),
     '3.0': array([1, 0])}
    >>> allIndices
    {'1.0': array([[0, 1]]),
     '2.0': array([[2, 2]]),
     '3.0': array([[1, 0], [1, 1], [1, 2], [2, 0], [2, 1]])}
    >>> outletINFO
    {"1.0": ["nan", 0],
     "2.0": ["3.0", 3.0],
     "3.0": ["1.0", 1.0]}
    """
    # Create 1D map version of the distributed maps, in order to get indices right
    catchment1D = flat(catchmentMap)
    outletMap[np.isnan(catchmentMap)] = np.nan
    outlet1D    = flat(outletMap)

    # Correct the distance map for a few cells that are situated outside the
    # main catchment (ID = 1), so the shape matches all other files
    dist1D = np.zeros(catchment1D.shape)
    dist1D[dist1D == 0] = np.nan
    dist1D[catchment1D >= 1] = flat(distMap)

    # Convert outlet to float and remove zeros (if not already done) and extract
    # all outletIDs
    outletMap = outletMap.astype(float)
    outletMap[outletMap == 0] = np.nan
    allOutlets = np.unique(outletMap[~np.isnan(outletMap)])
    # Create empty dictionary to store data in
    outletIndices = {}
    allIndices = {}
    outletINFO = {}

    outletIndices_1D = {}
    allIndices_1D = {}

    # Loop through all outlets
    for outlet in allOutlets:
        # Create empty vector to store data in
        # Find the index corresponding to that outlet
        outletIndices[str(outlet)] = np.transpose(np.where(outletMap == outlet))[0]
        outletIndices_1D[str(outlet)] = np.transpose(np.where(outlet1D == outlet))[0]
        # Find the indices related to the subcatchment of that outlet (same ID)
        allIndices[str(outlet)] = np.transpose(np.where(catchmentMap == outlet))
        allIndices_1D[str(outlet)] = np.transpose(np.where(catchment1D == outlet))

    for outlet in allOutlets:
        if outlet != 1.0:
            print(outlet)
            # Find the indices of the neighbours of the outlet
            nb_cells = get_neighbours(outletIndices[str(outlet)], shape = catchmentMap.shape)

            # Save only cell indexes situated outside the subcatchment
            nb_cells = nb_cells[np.where(catchmentMap[tuple(nb_cells.T)] != outlet)]
            # Find index with the lowest elevation
            lowest = nb_cells[np.where(elevMap[tuple(nb_cells.T)] == min(elevMap[tuple(nb_cells.T)]))]
            # Get the catchmentID that is related to this index
            nb_catch = catchmentMap[tuple(lowest.T)].astype(str)
            # Extract the cell ID
            nb_catch = nb_catch[0]

            # Find the distance of the outlet of interest
            upper_ind = outletIndices[str(outlet)]
            upper_dist = distMap[upper_ind[0],upper_ind[1]]
            # Find the distance of the outlet of the neighbouring catchment
            lower_ind = outletIndices[nb_catch]
            lower_dist = distMap[lower_ind[0],lower_ind[1]]

            outletINFO[str(outlet)] = [nb_catch, upper_dist - lower_dist]
        else:
            outletINFO[str(outlet)] = ["nan", 0]

    return outletIndices_1D, allIndices_1D, outletINFO

if __name__ == "__main__":
    import pickle

    distMap      = readMap("map_distance.asc")
    catchmentMap = readMap("map_catchment.asc")
    outletMap    = readMap("map_outlet.asc")
    elevMap      = readMap("map_dem.asc")

    distMap      = correct_NA(distMap,      -9999)
    catchmentMap = correct_NA(catchmentMap, -9999)
    outletMap    = correct_NA(outletMap,    -9999)
    elevMap      = correct_NA(elevMap,      -9999)

    outletLoc, catchLoc, outletInfo = subcatch_to_outlet(distMap, catchmentMap, outletMap, elevMap)

    def save_object(obj, filename):
        with open(filename, 'wb') as output:
            pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

    save_object([outletLoc, catchLoc, outletInfo], 'routingInfo.pkl')

    # Usage
    # outletLoc, catchLoc, outletInfo = pickle.load( open( "dS2_routingInfo.pkl", "rb" ) ) 


