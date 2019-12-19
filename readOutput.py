import glob
import numpy as np
import xarray as xr

def write_memmap(fileName, data, c_start, c_end):
    """ Write chunk of output data to a memmap file. The indexes of the chunk
    are passed onto the filename for easy access later. Data values are stored
    as numpy.float32 values."""
    fileName = "{}_{}_{}.dat".format(fileName, c_start, c_end)
    fp = np.memmap(fileName, dtype='float32', mode='w+', shape=data.shape)
    fp[:] = data[:]
    del fp

def memmap_to_nc(self, var):
    """ Read all the generated memmap files and write those to a NetCDF file,
    where the latitude and longitude are extracted from the input maps. Please
    note that this function removes the memmap files."""

    # List all files corresponding to the required pattern
    files = glob.glob("{}/{}*.dat".format(self.outdir, var))

    # Sort the files by 
    first_index = [int(f.split('\\')[-1].split("_")[1]) for f in files]
    files = [x for _,x in sorted(zip(first_index,files))]

    # Create empty matrix with the correct shape and datatype
    c_data = np.zeros(shape=self.shape, dtype="float32")

    # Offset in starting point if the start of sim period did not match the
    # start of the data period
    offset = int(files[0].split('\\')[-1].split("_")[1])

    # Loop through all files and import information into memory
    for file in files:
        info  = file.split('\\')[-1].split("_")
        start = int(info[1]) - offset
        stop   = int(info[2].replace(".dat","")) - offset
        shape = (stop-start, self.shape[1])
        data = np.memmap(file, dtype="float32", mode="r", shape = shape)

        c_data[start:stop] = data
        del data

    # Create empty matrix with correct 3D dimensions, based on the catchment
    # map and replace all values with NaNs
    allData = np.zeros(shape = (len(c_data), self.catchment.shape[0],
                                self.catchment.shape[1]), dtype = "float32")
    allData *= np.nan
    # For each timestep, put the data in the required shape
    for i in range(len(c_data)):
        allData[i][~np.isnan(self.catchment)] = c_data[i]
    # Transform 3D array into a dataset with correct dimensions
    # TODO: Fix the x and y positioning of the NetCDF files (now the indices are used)
    ncData = xr.Dataset({"val": (["time", "x", "y"],  allData)},
                        coords = {"time": self.sim_period,
                                  "x" : range(self.catchment.shape[0]),
                                  "y" : range(self.catchment.shape[1])})
    # Write to NetCDF file
    ncData.to_netcdf("{}/{}.nc".format(self.outdir, var))