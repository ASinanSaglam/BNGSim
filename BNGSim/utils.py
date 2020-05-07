import os, subprocess
from distutils import spawn

def find_BNG_path(BNGPATH=None):
    # Let's keep up the idea we pull this path from the environment
    if BNGPATH is None:
        try:
            BNGPATH = os.environ["BNGPATH"]
        except:
            pass
    # if still none, try pulling it from cmd line
    if BNGPATH is None:
        bngexec = "BNG2.pl"
        if test_bngexec(bngexec):
            print("BNG2.pl seems to be working")
            # get the source of BNG2.pl
            BNGPATH = spawn.find_executable("BNG2.pl")
            BNGPATH, _ = os.path.split(BNGPATH)
    else:
        bngexec = os.path.join(BNGPATH, "BNG2.pl")
        if test_bngexec(bngexec):
            print("BNG2.pl seems to be working")
        else:
            print("BNG2.pl not working, simulator won't run")
    return BNGPATH, bngexec

def test_bngexec(bngexec):
    rc = subprocess.run(["perl",bngexec])
    if rc.returncode == 0:
        return True
    else:
        return False

