import os
import glob
# import time
import numpy as np
import pandas as pd
from tqdm import trange

import readInput as rI
import funcFluxes as fF
import funcSolver as fS
import funcRouting as fR
import funcParams as fP
import readOutput as rO

class dS2:

    def model_setup_class(self, input_class):
        ''' Function to read the settings class, read the required input files,
        set the parameter values and set the model options. '''
        print(">>> Preparing input...")
        import copy as cp
        settings = cp.deepcopy(input_class)
        rI.read_class(self, settings)

    def change_param(self, par_name, value):
        ''' Use this function to change parameter values, as some parameters
        are corrected to different units to be used in the model. This function
        ensures that the parameter will have the correct units after changing
        the value. '''
        fP.set_parameter(self, par_name, value)

    def __read_chunk(self, chunk_end, offset, i):
        ''' Function to read part of the forcing data (in memmap format) and
        calculate all other forcing related variables, to reduce computational
        time. '''
        chunk_start = chunk_end
        chunk_end = min((self.tsteps + offset) - chunk_start, self.chunk_size) + i + offset

        # Read precipitation data
        self.P = rI.read_memmap(self, fileName=self.P_file)[chunk_start:chunk_end]
        self.ET_CHUNK = rI.read_memmap(self, fileName=self.ET_file)[chunk_start:chunk_end]
        # Read chunks of the required input data for snow, set Rs to zero if not present
        if self.SnowFLAG:
            T  = rI.read_memmap(self, fileName=self.T_file) [chunk_start:chunk_end]
            try:
                Rs = rI.read_memmap(self, fileName=self.Rs_file)[chunk_start:chunk_end]
            except FileNotFoundError:
                Rs = 0
            self.Prain_CHUNK = np.where(T > self.T0, self.P, 0)
            self.Psnow_CHUNK = np.where(T <= self.T0, self.P, 0)

            # Reshape parameter
            if type(self.ddf) == np.ndarray:
                ddf = np.tile(self.ddf, len(T)).reshape(T.shape)
            else:
                ddf = self.ddf

            self.Smelt_POT = fF.snowmelt_potential(T = T, Rg = Rs, T0 = self.T0,
                                              ddf = ddf, fact = self.rdf)
        else:
            self.Prain_CHUNK = self.P
            self.Psnow_CHUNK = np.zeros(self.P.shape, dtype=self.dtype)
        return chunk_end

    def __write_tmp_output(self, var, chunk_start, chunk_end):
        ''' Write temporary output files, as result of the chunking
        functionality. This will only be done for the variables that are
        demanded as output. '''
        if var == "WatBal":
            P_sum = np.mean(self.P)
            ET_sum = np.mean(self.ETpot) * self.eps
            Qsim_sum = np.mean(self.Qsim)
            #TODO: WatBal does not include storage
            self.WatBal = self.WatBal.append(pd.Series([P_sum, ET_sum, Qsim_sum],
                                                       index=self.WatBal.columns),
                                            ignore_index=True)
        else:
            fileName = "{}/{}".format(self.outdir, var)
            dat = getattr(self, var)
            rO.write_memmap(fileName, dat, chunk_start, chunk_end)

    def generate_runoff(self, progress=False):
        '''Loops through all the time steps to calculate the discharge. '''
        print(">>> Generating discharge...")

        # reshape the parameters (if only a single value), to allow for slicing
        if self.evap_reduction:
            self.reduction_idx = []
            for par in ["alpha", "beta", "gamma"]:
                tmp_par = getattr(self, par)

                # test whether parameters are an array or single value
                if type(tmp_par) != np.ndarray:
                    # if not an array, reshape the parameter to the correct shape
                    tmp_par = np.repeat(tmp_par, self.shape[1])
                    setattr(self, par, tmp_par)

        #======================================================================
        # Remove all contents in output directory
        #======================================================================
        try:
            for file_name in os.listdir(self.outdir):
                os.remove("{}/{}".format(self.outdir, file_name))
        except FileNotFoundError:
            os.mkdir(self.outdir)

        #======================================================================
        # Find forcing index corresponding to starting date
        #======================================================================
        offset = int(np.where(self.data_period==self.sim_period[0])[0])
        if offset + len(self.sim_period) > len(self.data_period):
            raise ValueError("Simulation period is too long")

        #======================================================================
        # Read first input chunk
        #======================================================================
        chunk_start = 0 + offset
        chunk_end = self.__read_chunk(chunk_end = chunk_start,
                                      offset = offset, i = 0)


        #======================================================================
        # Prepare output files
        #======================================================================
        allVars = ["WatBal", "Qsim", "Smelt", "Sstore"]
        outputVars = []
        for var in allVars:
            if getattr(self, "{}_flag".format(var)):
                if (var == "Smelt" or var == "Sstore" and self.SnowFLAG):
                    outputVars += [var]
                elif (var != "Smelt" and var != "Sstore"):
                    outputVars += [var]
                if var != "WatBal":
                    setattr(self, var, np.zeros(self.P.shape, dtype = self.dtype))
                else:
                    self.WatBal = pd.DataFrame(columns = ("P", "ET", "Qsim"))
                    setattr(self, "ETpot", np.zeros(self.P.shape, dtype = self.dtype))
        self.outputVars = outputVars

        #======================================================================
        # Initialize values
        #======================================================================
        # Set initial discharge
        Qsim         = np.broadcast_to(self.init_Qsim, (self.shape[1],))
        # Set initial snow storage
        if self.SnowFLAG:
            Sstore = self.init_Sstore

        # Set the for loop to be a tqdm loop if progress bar is required
        if progress == True:
            sim_length = trange(0, self.tsteps)
        else:
            sim_length = range(0, self.tsteps)

        # Chunk indexing variable, set to -1 so it is zero in the first loop
        j = -1
        # Variable to store the number of internal timestep per timestep
        self._num_dt = []
        # self.solver_time = []

        for i in sim_length:
            # Set index for the forcing data (since it is chunked)
            j += 1
            #==================================================================
            # Output the results if required
            #==================================================================
            if i%self.chunk_size == 0:
                if i != 0:
                    for var in outputVars:
                        self.__write_tmp_output(var, chunk_start, chunk_end)

                    #==========================================================
                    # Read new chunk of forcing data
                    #==========================================================
                    j = 0
                    chunk_start = chunk_end
                    chunk_end = self.__read_chunk(chunk_end = chunk_end, offset = offset, i = i)
                    #==========================================================
                    # Create new empty arrays for the output variables
                    #==========================================================
                    for var in outputVars:
                        setattr(self, var, np.zeros(self.P.shape, dtype = self.dtype))

            #==================================================================
            # Run the model
            #==================================================================
            # Solid and liquid precipitation
            Prain = self.Prain_CHUNK[j]
            Psnow = self.Psnow_CHUNK[j]

            # Evaporation
            ETpot = self.ET_CHUNK[j]

            # Calculate the snowmelt and update snow storage
            if self.SnowFLAG:
                Smelt_mmdt = np.minimum(self.Smelt_POT[j] * self.dt, Sstore)
                Sstore = Sstore + self.dt * Psnow - Smelt_mmdt
                Smelt = Smelt_mmdt / self.dt
            else:
                Smelt = 0

            # Calculate the discharge generated in each cell
            # start = time.time()
            Qnew, num_dt = fS.flexible_RK4(P = Prain + Smelt, ET = ETpot * self.eps, Q = Qsim,
                                alpha = self.alpha, beta = self.beta, gamma = self.gamma,
                                dt=self.dt, LB=self.LB,
                                gQ_mdiff=2.0, dt_reduction=.15, min_dt=10, max_dt=50)
            # Qnew, num_dt = fS.cashkarp_wrapper(P = Prain + Smelt, ET = ETpot * self.eps, Q = Qsim,
            #                 alpha = self.alpha, beta = self.beta, gamma = self.gamma,
            #                 dt = 1, LB = 1e-10, dt_factor = 1000, min_dt = 5, max_dt = 100)
            # stop = time.time()
            # self.solver_time.append(stop-start)

            # Perform evaporation reduction (set ET = 0 if Q < Qt)
            if self.evap_reduction == True:
                Qnew, idx = fF.evaporation_reduction(QwithET=Qnew, Q_init=Qsim, P=Prain + Smelt,
                                                     alpha=self.alpha, beta=self.beta, gamma=self.gamma,
                                                     Q_threshold=self.Qt, LB=self.LB, dt = self.dt)
                self.reduction_idx.append(idx)

            Qsim = Qnew
            self._num_dt.append(num_dt)

            #==================================================================
            # Set the value to a class variable
            #==================================================================
            for var in outputVars:
                if var == "WatBal":
                    getattr(self, "ETpot")[j] = locals()["ETpot"]
                else:
                    getattr(self, var)[j] = locals()[var]
        #======================================================================
        # Save the output variables to a file
        #======================================================================
        for var in outputVars:
            self.__write_tmp_output(var, chunk_start, chunk_end)

        #======================================================================
        # Create correct water balance file
        #======================================================================
        if self.WatBal_flag:
            self.WatBal.loc["total"] = self.WatBal.sum()
            diff = self.WatBal.loc["total","P"] - self.WatBal.loc["total","ET"] - self.WatBal.loc["total","Qsim"]
            self.WatBal.loc["total","diff"] = diff
            self.WatBal.to_csv("{}/WaterBalance.csv".format(self.outdir))


    def export(self):
        '''Export the output variables from numpy memmap files to a NetCDF, and
        deletes the numpy memmap files (except for the Qsim variable, as this
        data is required for the routing function. '''
        print(">>> Writing output...")
        for var in self.outputVars:
            if var == "Qsim":
                if getattr(self, "{}_flag_nc".format(var).format(var)):
                    rO.memmap_to_nc(self, var)
            else:
                rO.memmap_to_nc(self, var)
                # Remove the .dat files
                for f in glob.glob("{}/{}*.dat".format(self.outdir, var)) :
                    os.remove(f)


    def routing_with_diffusion(self, main_ID=1.0, delete=True, trackwater=False):
        print(">>> Routing with diffusion...")

        # Create a list of all outlets
        allOutlets = list(map(float, self.outletLoc.keys()))
        # Create empty dictionary to store results in
        self.Qrout = {}
        for outlet in allOutlets:
            self.Qrout[str(outlet)] = np.repeat(0., self.shape[0] +
                                      int(np.nanmax(self.dist1D) * self.tau))

        if trackwater:
            # Create empty dataframe to store the routed water per timestep per basin in
            # self.Qorigin = pd.DataFrame(0, columns=range(self.shape[1]),
            #                                index=np.repeat(0., self.shape[0]))
            self.Qorigin = pd.DataFrame(0, columns=np.array(allOutlets).astype(str),
                                           index=range(self.shape[0]+int(np.nanmax(self.dist1D) * self.tau)))
            # self.track_count = self.Qorigin.copy()

        # Count number of pixels contributing to each basin
        count = {str(x): len(self.catchLoc[str(x)]) for x in allOutlets}
        for outlet in self.outletInfo:
            dwn_outlet = self.outletInfo[outlet][0]
            while dwn_outlet != "nan":
                count[dwn_outlet] += len(self.catchLoc[outlet])
                dwn_outlet = self.outletInfo[dwn_outlet][0]

        # Define the maximum window size
        max_wid = int(np.max(self.dist1D) * self.tau * self.lag_to_window)
        max_lag = int(np.max(self.dist1D) * self.tau)
        max_window =  max_lag + max_wid

        if max_lag < max_wid//2:
            slice_offset = max_wid//2 - max_lag
        else:
            slice_offset = 0

        # Find the filenames of the memmap files and sort
        files = glob.glob("{}/{}*.dat".format(self.outdir, "Qsim"))
        first_index = [int(f.split('\\')[-1].split("_")[1]) for f in files]
        files = [x for _,x in sorted(zip(first_index,files))]
        # Find the offset of the data
        offset = int(files[0].split('\\')[-1].split("_")[1])

        copy_end = {}
        # Loop through all files and import information into memory
        for file in files:
            # Extract chunk information and set the indices
            info  = file.split('\\')[-1].split("_")
            start = int(info[1]) - offset
            stop   = int(info[2].replace(".dat","")) - offset
            shape = (stop-start, self.shape[1])
            # Read the memmap data
            data = np.memmap(file, dtype=self.dtype, mode="r", shape = shape)

            if trackwater:
                # Create temporary dataframe, used in mutliOutlet_routing()
                self.Qorigin_tmp = pd.DataFrame(0,
                    columns=np.array(allOutlets).astype(str),
                    index=range(shape[0] + int(np.nanmax(self.dist1D) * self.tau * 3)))
                # self.track_count_tmp = self.Qorigin_tmp.copy()
                # self.Qorigin_tmp = pd.DataFrame(0, columns=range(self.shape[1]),
                #                            index=range(len(np.repeat(0., shape[0] +
                #                         int(np.nanmax(self.dist1D) * self.tau * 2)))))

            if file == files[0]:
                # Rout the chunk of discharge data
                Qtmp = fR.with_diffusion(self, Qsim = data, main_ID=1.0, trackwater=trackwater)
                last = np.array(data)[-max_window:]
            else:
                data = np.insert(data, 0, last, axis=0)
                Qtmp = fR.with_diffusion(self, Qsim = data, main_ID=1.0, trackwater=trackwater)
                last = np.array(data)[-max_window:]

            # Write temporary values to the total dictionary
            for outlet in Qtmp:
                outlet = str(outlet)
                # Slice the trailing zeroes, to prevent replacing errors
                val = np.trim_zeros(Qtmp[outlet], trim="b")
                if file != files[0]:
                    val = val[max_window-slice_offset:]
                    val = val[:self.chunk_size]
                    self.Qrout[outlet][copy_end[outlet]-slice_offset:copy_end[outlet]-slice_offset + len(val)] = val/count[outlet]

                    if trackwater:
                        tmp_val = np.trim_zeros(self.Qorigin_tmp.loc[:, outlet].values, trim="b")
                        tmp_val = tmp_val[max_window-slice_offset:]
                        tmp_val = tmp_val[:self.chunk_size-1]
                        self.Qorigin.loc[copy_end[outlet]-slice_offset:copy_end[outlet]-slice_offset+len(tmp_val)-1,outlet] = tmp_val

                    copy_end[outlet] = copy_end[outlet] + len(val)
                else:

                    if trackwater:
                        # tmp_val = self.Qorigin_tmp.loc[:self.chunk_size-1,outlet]
                        tmp_val = np.trim_zeros(self.Qorigin_tmp.loc[:,outlet], trim="b")
                        tmp_val = tmp_val[:self.chunk_size]
                        self.Qorigin.loc[:len(tmp_val)-1, outlet] = tmp_val

                    val = val[:self.chunk_size]
                    self.Qrout[outlet][:len(val)] = val/count[outlet]
                    copy_end[outlet] = len(val)
            # # Add values to the complete Qorigin dataframe and remove the temporary file
            # if trackwater:
            #     # self.Qorigin[start:start+len(self.Qorigin_tmp)] += \
            #     self.Qorigin.iloc[start:start + len(self.Qorigin_tmp)] += \
            #         self.Qorigin_tmp.iloc[:self.tsteps - start].values
            #         # self.Qorigin_tmp.iloc[:min(self.tsteps, start+len(self.Qorigin_tmp))].values

            #     self.track_count.iloc[start:start + len(self.Qorigin_tmp)] += \
            #         self.track_count_tmp.iloc[:self.tsteps - start].values
            #         # self.track_count_tmp.iloc[:min(self.tsteps, start+len(self.track_count_tmp))].values
            #     # del self.Qorigin_tmp, self.track_count_tmp

            # Close and delete the memmap file
            del data#TODO, self.Qorigin_tmp
            if delete:
                os.remove(file)

        # Correct the dictionary to the total simulation length
        for key in self.Qrout.keys():
            self.Qrout[key] = self.Qrout[key][:self.tsteps]

        if trackwater:
            self.Qorigin = self.Qorigin[:self.tsteps] / self.shape[1]

    def routing(self, main_ID=1.0, delete=True, trackwater=False):
        '''Transports water from each pixel to the outlet(s), applying a time
        lag based on the distance of each pixel to the outlet(s). '''
        print(">>> Routing...")

        # Create a list of all outlets
        allOutlets = list(map(float, self.outletLoc.keys()))
        # Create empty dictionary to store results in
        self.Qrout = {}
        for outlet in allOutlets:
            self.Qrout[str(outlet)] = np.repeat(0., self.shape[0] +
                                      int(np.nanmax(self.dist1D) * self.tau))
        if trackwater:
            # Create empty dataframe to store the routed water per timestep per basin in
            self.Qorigin = pd.DataFrame(0,
                                        columns=np.array(allOutlets).astype(str),
                                        index=range(len(np.repeat(0., self.shape[0] +
                                                    int(np.nanmax(self.dist1D) * self.tau * 2)))))

        # Find the filenames of the memmap files and sort
        files = glob.glob("{}/{}*.dat".format(self.outdir, "Qsim"))
        first_index = [int(f.split('\\')[-1].split("_")[1]) for f in files]
        files = [x for _,x in sorted(zip(first_index,files))]
        # Find the offset of the data
        offset = int(files[0].split('\\')[-1].split("_")[1])

        # Loop through all files and import information into memory
        for file in files:
            # Extract chunk information and set the indices
            info  = file.split('\\')[-1].split("_")
            start = int(info[1]) - offset
            stop   = int(info[2].replace(".dat","")) - offset
            shape = (stop-start, self.shape[1])
            # Read the memmap data
            data = np.memmap(file, dtype=self.dtype, mode="r", shape = shape)

            if trackwater:
                # Create temporary dataframe, used in mutliOutlet_routing()
                self.Qorigin_tmp = pd.DataFrame(0,
                    columns=np.array(allOutlets).astype(str),
                                        index=range(len(np.repeat(0., shape[0] +
                                        int(np.nanmax(self.dist1D) * self.tau * 2)))))

            # Rout the chunk of discharge data
            Qtmp = fR.multiOutlet_routing(self, data, main_ID, trackwater)

            # Write temporary values to the total dictionary
            for outlet in Qtmp:
                # Slice the trailing zeroes, to prevent replacing errors
                val = np.trim_zeros(Qtmp[str(outlet)], trim="b")
                self.Qrout[str(outlet)][start:start+len(val)] += val

            # Add values to the complete Qorigin dataframe and remove the temporary file
            if trackwater:
                self.Qorigin[start:start+len(self.Qorigin_tmp)] += \
                    self.Qorigin_tmp[:min(self.tsteps, start+len(self.Qorigin_tmp))].values
                del self.Qorigin_tmp

            # Close and delete the memmap file
            del data
            if delete:
                os.remove(file)

        # Correct the dictionary to the total simulation length
        for key in self.Qrout.keys():
            self.Qrout[key] = self.Qrout[key][:self.tsteps]

        if trackwater:
            # Slice to correct length, and correct for the number of pixels
            self.Qorigin = self.Qorigin[:self.tsteps]
            self.Qorigin /= len(self.catchment[~np.isnan(self.catchment)])
            self.Qorigin.index = pd.DatetimeIndex(self.sim_period)



    def _routing_corrected_mm(self, fname_size_map, main_ID=1.0):
        ''' Same as normal routing function, but correts for the size of the
        pixels. USE ONLY WHEN PIXELS DO NOT HAVE A UNIFORM SIZE'''
        print(">>> Routing...")

        # Create a list of all outlets
        allOutlets = list(map(float, self.outletLoc.keys()))
        # Create empty dictionary to store results in
        self.Qrout = {}
        for outlet in allOutlets:
            self.Qrout[str(outlet)] = np.repeat(0., self.shape[0] +
                                      int(np.nanmax(self.dist1D) * self.tau))

        # Find the filenames of the memmap files and sort
        files = glob.glob("{}/{}*.dat".format(self.outdir, "Qsim"))
        first_index = [int(f.split('\\')[-1].split("_")[1]) for f in files]
        files = [x for _,x in sorted(zip(first_index,files))]
        # Find the offset of the data
        offset = int(files[0].split('\\')[-1].split("_")[1])

        ##### Read the pixel size map #####
        self.size = np.loadtxt("{}/{}".format(self.indir, fname_size_map), skiprows=0)
        self.size[self.size==-9999] = np.nan
        self.size = self.size[~np.isnan(self.size)]
        # total_size = sum(self.size)
        mean_size = np.mean(self.size)

        # Loop through all files and import information into memory
        for file in files:
            # Extract chunk information and set the indices
            info  = file.split('\\')[-1].split("_")
            start = int(info[1]) - offset
            stop   = int(info[2].replace(".dat","")) - offset
            shape = (stop-start, self.shape[1])
            # Read the memmap data
            data = np.memmap(file, dtype=self.dtype, mode="r", shape = shape)
            # Correct data for the pixel size
            data = (data * self.size) / mean_size

            # Rout the chunk of discharge data
            Qtmp = fR.multiOutlet_routing(self, data, main_ID)

            # Write temporary values to the total dictionary
            for outlet in Qtmp:
                # Slice the trailing zeroes, to prevent replacing errors
                val = np.trim_zeros(Qtmp[str(outlet)], trim="b")
                self.Qrout[str(outlet)][start:start+len(val)] += val

            # Close and delete the memmap file
            del data
            os.remove(file)
