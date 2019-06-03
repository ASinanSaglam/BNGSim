import os, subprocess

def find_BNG_path(BNGPATH):
    # Let's keep up the idea we pull this path from the environment
    if BNGPATH == "":
        try:
            BNGPATH = os.environ['BNGPATH']
        except KeyError:
            print("BNGPATH is not given or set in the environment")
    # Raise a warning if we don't have access to BNGPath now
    if BNGPATH == "":
        bngexec = "BNG2.pl"
    else:
        bngexec = BNGPATH + "/BNG2.pl"
    if not test_bngexec(bngexec):
        print("BNG2.pl not working, simulator won't run")
    else:
        print("BNG2.pl seems to be working")
    return BNGPATH, bngexec

def test_bngexec(bngexec):
    rc = subprocess.run([bngexec])
    if rc.returncode == 0:
        return True
    else:
        return False

