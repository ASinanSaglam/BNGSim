import os, sys, subprocess, h5py, pickle
from multiprocessing import Pool
from shutil import copyfile
import numpy as np

from BNGSimulator import BNGSimulator
from BNGResult import BNGResult

class BNGSim:
    """
    The starting goals of this class is as follows:
    
    - Takes in a path, either existing or not, and associates itself with that path
    - Can be given a BNGL file and BNGPATH to the BioNetGen installation. If so
    this class then can run the said BNG simulation (which will be delagated to 
    another class, this one just sets it all up)
    - This can be skipped, in which case it should assume this is for analysis. Then
    it automatically parses for existing files and loads them in. 
    
    TODO: Actually manipulate the BNGL file so we can setup the actions block if nothing else
    TODO: Search for a net file, if exists fall back to run_network/nfsim binaries and don't call
    BNG2.pl (what other file do we need for nfsim again?). This essentially makes it 
    usable by WESTPA directly.
    TODO: Incorporate PySB simulator over the super basic simulator we have here. Similar for 
    BNGL manipulation too.
    """
    def __init__(self, path, BNGPATH="", bngl=None, ncores=1, cleanup=True):
        # Basic path setting here
        self._setup_working_path(path)
        self.ncores = ncores
        self.bngl = bngl
        self.BNGPATH = BNGPATH
        # Initialize results
        self.results = []
        self.combined_results = None
        self.cleanup = cleanup

    def _setup_simulator(self, bngl, BNGPATH):
        # Simulator WILL assume we are in the correct folder so make sure we are by this point
        simulator = BNGSimulator(bngl, BNGPATH)
        self.bngl = simulator.bngl
        self.BNGPATH = simulator.BNGPATH
        self.bngexec = simulator.bngexec
        return simulator
    
    def _setup_simulators(self, bngl, BNGPATH, ncores):
        paths = [os.path.join(self.path, "{:08d}".format(i)) for i in range(ncores)]
        simulators = [BNGSimulator(bngl, BNGPATH, path, cleanup) for path in paths]
        self.bngls = [simulator.bngl for simulator in simulators]
        self.BNGPATH = simulators[0].BNGPATH
        self.bngexec = simulators[0].bngexec
        return simulators
    
    def _call_into_simulator(self, simulator):
        res = simulator.run()
        return res

    def run_simulation(self, nsims=1):
        # Setting stuff up for simulation in case we want that
        if nsims == 1:
            self.simulator = self._setup_simulator(self.bngl, self.BNGPATH)
        else:
            self.simulators = self._setup_simulators(self.bngl, self.BNGPATH, nsims)

        self.ensure_working_path()
        if self.ncores == 1:
            if not self.simulator.is_ready():
                print("Can't run a simulation because the simulator is not setup correctly")
                return
            else:
                for i in range(nsims):
                    result = self.simulator.run()
                    if result is not None:
                        result.set_name("simulation_{:08d}".format(len(self.results)))
                        self.results.append(result)
        elif self.ncores > nsims:
            print("running parallel with {} cores".format(self.ncores))
            p = Pool(nsims)
            para_res = p.map(self._call_into_simulator, self.simulators)
            for res in para_res:
                res.set_name("simulation_{:08d}".format(len(self.results)))
                self.results.append(res)
        else:
            print("running parallel with {} cores".format(self.ncores))
            p = Pool(self.ncores)
            para_res = p.map(self._call_into_simulator, self.simulators)
            for res in para_res:
                res.set_name("simulation_{:08d}".format(len(self.results)))
                self.results.append(res)
    
    def _setup_working_path(self, path):    
        if not os.path.isdir(path): 
            print("Given simulation path does not exist or invalid, trying to create")
            try: 
                os.mkdir(path)
            except OSError:
                print("Failed to create given path")
                os.exit(1)
        
        self.path = path
        # Switch there 
        if os.getcwd() != self.path:
            os.chdir(self.path)
        return     
   
    def ensure_working_path(self, path=None):
        if os.getcwd() != self.path:
            os.chdir(self.path)
        return 
    
    def collect_results(self):
        self.results.append(BNGResult(os.getcwd(),""))
        return
    
    def save_results(self, fname="results.h5", combined=False):
        """
        Saves results in an hdf5 file
        
        TODO: How do we want to tackle combined results? Do we want this to 
        automatically save that instead? Or do we want to have a flag for them?
        """
        dt = h5py.special_dtype(vlen=str)
        with h5py.File(fname, "w") as h:
            for result in self.results:
                grp = h.create_group(result.name)
                ## Look up how to save strings here
                #grp.create_dataset("bngl_file", data=result.bngl, dtype="S10")
                #if isinstance(result.bngl, dict):
                #    bngl_to_save = result.bngl[result.bngl.keys[0]]
                #    grp.attrs["bngl_file"] = bngl_to_save
                #else:
                #    grp.attrs["bngl_file"] = result.bngl

                #if isinstance(result.net, dict):
                #    net_to_save = result.bngl[result.bngl.keys[0]]
                #    grp.attrs["net_file"] = net_to_save
                #else:
                #    grp.attrs["net_file"] = result.net
                # Let's save cdat/gdats
                if isinstance(result.gdat, dict):
                    gdat_grp = grp.create_group("gdat")
                    for key in result.gdat.keys():
                        gdat_grp.create_dataset(key, data=result.gdat[key], dtype=result.gdat[key].dtype)
                else:
                    gdat_obj = grp.create_dataset("gdat", data=result.gdat, dtype=result.gdat.dtype)
                # now the same for cdat
                if isinstance(result.cdat, dict):
                    cdat_grp = grp.create_group("cdat")
                    for key in result.cdat.keys():
                        cdat_grp.create_dataset(key, data=result.cdat[key], dtype=result.cdat[key].dtype)
                else:
                    cdat_obj = grp.create_dataset("cdat", data=result.cdat, dtype=result.cdat.dtype)
            if combined:
                if self.combined_results is not None:
                    combined_obj = h.create_dataset("combined_results", data=self.combined_results,
                                                   dtype=self.combined_results.dtype)
    
    def combine_results(self):
        """
        Combines all results gdat arrays into a single array. Ensure that the dtype for all
        results is the same! (so the same set of observables but can be different lengths)
        
        TODO: How do we really want to tackle this? Only gdat? Include cdat? Try to combine
        everything? 
        """
        if len(self.results) == 0:
            print("combine_results is called without any results loaded in")
            return None
        
        # We'll use the same dtype and use the maximum length gdat array
        # if all the same, fine, if not the value will be set to NaN
        # The DTYPE has to be the same for all for this to work!
        nres = len(self.results)
        max_len = max([self.results[i].gdat.shape[0] for i in range(nres)])
        self.combined_results = np.empty((nres,max_len), dtype=self.results[0].gdat.dtype)
        self.combined_results[:] = np.nan
        for i, result in enumerate(self.results):
            self.combined_results[i] = result.gdat[:]
        print("Results combined in combined_results attribute")