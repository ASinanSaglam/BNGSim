import os, sys, subprocess, h5py, pickle
from multiprocessing import Pool
from shutil import copyfile
from shutil import rmtree

import numpy as np

from BNGResult import BNGResult

class BNGSimulator:
    """
    Assumes it's in the correct directory, takes the BNGL file (or eventually
    some other structure based on BNGL) and copies into the directory if necessary
    and runs BNG2.pl. This is the basic method of operation at the moment, later will
    either be completely replaced by PySB simulator object or other stuff added
    """
    def __init__(self, bngl, BNGPATH, path=None, cleanup=True):
        # If path is given we assume we are doing a multiprocess run
        if path is not None:
            self._setup_working_path(path)
        self.set_BNG_path(BNGPATH)
        self.setup_BNGL(bngl)
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
        return

    def ensure_working_path(self):
        print("Ensuring path: {}".format(self.path))
        if os.getcwd() != self.path:
            os.chdir(self.path)

    def is_ready(self):
        """
        Checks if the simulator is ready to run
        """
        if not self.test_bngexec():
            return False
        if not os.path.isfile(self.bngl):
            return False
        return True
    
    def setup_BNGL(self, bngl):
        if bngl is None:
            # Let's see if we can divine the BNGL file first
            bngl_files = list(filter(lambda x: x.endswith(".bngl"), os.listdir(os.getcwd())))
            if len(bngl_files) == 0:
                print("No BNGL file given or exists in the PWD, simulator won't run")
            elif len(bngl_files) > 1:
                print("More than one BNGL file exists in PWD, using {}".format(bngl_files[0]))
            else:
                print("BNGL file found in PWD: {}".format(bngl_files[0]))
                self.bngl = bngl_files[0]
        else:
            if not os.path.isfile(bngl):
                print("BNGL file given doesn't exist, simulator won't run")
            else:
                # Now assuming BNGL file given exists, we are in the right folder
                # self.bngl will point to the copied 
                self.bngl_name = os.path.basename(bngl)
                self.bngl = os.path.join(os.getcwd(), self.bngl_name)
                if not os.path.isfile(self.bngl):
                    try:
                        copyfile(bngl, self.bngl)
                    except IOError:
                        print("BNGL file can't be copied to PWD, simulator won't run")
        return 
    
    def set_BNG_path(self, BNGPATH):
        # Let's keep up the idea we pull this path from the environment
        if BNGPATH == "":
            try:
                BNGPATH = os.environ['BNGPATH']
            except KeyError:
                print("BNGPATH is not given or set in the environment")
        # Raise a warning if we don't have access to BNGPath now
        self.BNGPATH = BNGPATH
        if BNGPATH == "":
            self.bngexec = "BNG2.pl"
        else:
            self.bngexec = self.BNGPATH + "/BNG2.pl"
        if not self.test_bngexec():
            print("BNG2.pl not working, simulator won't run")
        else:
            print("BNG2.pl seems to be working")
        return
    
    def test_bngexec(self):
        rc = subprocess.run([self.bngexec])
        if rc.returncode == 0:
            return True
        else:
            return False

    def clean_sim_folder(self):
        if os.getcwd() == self.path:
            rmtree(self.path)
            os.chdir("../")
        else:
            rmtree(self.path)
        print("Cleaned up simulation folder")
        
    def run(self):
        self.ensure_working_path()
        rc = subprocess.run([self.bngexec, self.bngl])
        #print(rc)
        if rc.returncode == 0:
            print("Simulation succesful, loading results")
            self.result = BNGResult(os.getcwd(), self.bngl)
            if self.cleanup:
                print("Cleanup asked, removing simulation folder {}".format(self.path))
                self.clean_sim_folder()
            return self.result
        else:
            print("Simulation failed")
            print(rc.stdout)
            print(rc.stderr)
            return None
