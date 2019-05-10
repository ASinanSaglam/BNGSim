import os, subprocess

def set_BNG_path(obj, BNGPATH):
    # Let's keep up the idea we pull this path from the environment
    if BNGPATH == "":
        try:
            BNGPATH = os.environ['BNGPATH']
        except KeyError:
            print("BNGPATH is not given or set in the environment")
    # Raise a warning if we don't have access to BNGPath now
    obj.BNGPATH = BNGPATH
    if BNGPATH == "":
        obj.bngexec = "BNG2.pl"
    else:
        obj.bngexec = obj.BNGPATH + "/BNG2.pl"
    if not test_bngexec(obj):
        print("BNG2.pl not working, simulator won't run")
    else:
        print("BNG2.pl seems to be working")
    return

def test_bngexec(obj):
    rc = subprocess.run([obj.bngexec])
    if rc.returncode == 0:
        return True
    else:
        return False

