import os, sys, subprocess, h5py, pickle
from multiprocessing import Pool
from shutil import copyfile
from shutil import rmtree

import numpy as np

import BNGUtils
from BNGResult import BNGResult

class BNGSimulator:
    """
    Assumes it's in the correct directory, takes the BNGL file (or eventually
    some other structure based on BNGL) and copies into the directory if necessary
    and runs BNG2.pl. This is the basic method of operation at the moment, later will
    either be completely replaced by PySB simulator object or other stuff added
    """
    def __init__(self, parser, BNGPATH, path=None, cleanup=True):
        # If path is given we assume we are doing a multiprocess run
        self.parser = parser
        if path is not None:
            self._setup_working_path(path)
        self.prep_bngl()
        BNGUtils.set_BNG_path(self, BNGPATH)
        self.result = None
        self.cleanup = cleanup
        return
    
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
        # Get our 
        return

    def prep_bngl(self):
        self.bngl_path = self.parser.export_bngl()

    def ensure_working_path(self):
        print("Ensuring path: {}".format(self.path))
        if os.getcwd() != self.path:
            os.chdir(self.path)

    def is_ready(self):
        """
        Checks if the simulator is ready to run
        """
        if not BNGUtils.test_bngexec(self):
            return False
        if not os.path.isfile(self.bngl_path):
            return False
        return True
    
    def clean_sim_folder(self):
        if os.getcwd() == self.path:
            rmtree(self.path)
            os.chdir("../")
        else:
            rmtree(self.path)
        print("Cleaned up simulation folder")
        
    def run(self):
        self.ensure_working_path()
        rc = subprocess.run([self.bngexec, self.bngl_path])
        #print(rc)
        if rc.returncode == 0:
            print("Simulation succesful, loading results")
            self.result = BNGResult(os.getcwd(), self.bngl_path)
            if self.cleanup:
                print("Cleanup asked, removing simulation folder {}".format(self.path))
                self.clean_sim_folder()
            return self.result
        else:
            print("Simulation failed")
            print(rc.stdout)
            print(rc.stderr)
            return None
