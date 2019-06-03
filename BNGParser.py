import os, subprocess, copy
import BNGUtils
from bs4 import BeautifulSoup as BS

# TODO: Collate together useful functions and make them into BNGUtils 
# instead of copying them into classes 

class BNGParser:
    '''
    This is a very simple BNGL reader that's able to read,y write the file, 
    strip actions that are in the file originally in there and able to append
    lines to the file. 
    '''
    def __init__(self, bngl, BNGPATH, run_params=None):
        if bngl is None:
            self.find_and_set_bngl()
        else:
            self.bngl_file = bngl
        if self.bngl_file is not None:
            if not run_params is None:
                self.init_bngl(self.bngl_file, no_action=True)
            else:
                self.init_bngl(self.bngl_file)
        if BNGPATH != "": 
            # We assume BNGPATH is legit and we'll use it for 
            # XML generation
            BNGPATH, bngexec = BNGUtils.find_BNG_path(BNGPATH)
            self.BNGPATH = BNGPATH
            self.bngexec = bngexec 
            self.gen_and_load_XML()
        if not run_params is None:
            self.run_params = run_params
            self.add_params_to_bngl()

    def gen_and_load_XML(self):
        self.clean_actions()
        rc = subprocess.run([self.bngexec, "--xml", self.bngl_file])
        bname = os.path.basename(self.bngl_file)
        mname = bname.replace(".bngl", "")
        xmlName = mname + ".xml"
        if rc.returncode == 0:
            print("XML generated")
            self.loaded_xml = self.load_xml(xmlName)
            print("XML loaded")

    def load_xml(self, xml):
        f = open(xml, 'r')
        sxml = BS(f, 'xml')
        f.close()
        return sxml

    def init_bngl(self, bngl_file, no_action=False):
        self.bngl_list = self.read_to_list(bngl_file)
        self.bngl = "".join(self.bngl_list)
        if no_action:
            self.clean_actions()
        self.export_bngl("current.bngl")
        # TODO: This bit is a bit hacky, gotta figure out a better way 
        self.bngl_file = "current.bngl"

    def read_to_list(self, bngl_file):
        with open(bngl_file, 'r') as f:
            bngl_list = f.readlines()
        return bngl_list

    def clean_actions(self):
        no_action_bngl = []
        modified = False
        for line in self.bngl_list:
            if "end model" in line: 
                no_action_bngl.append(line)
                modified = True
                break
            no_action_bngl.append(line)
        if modified:
            self.bngl_list = no_action_bngl
            self.original_bngl = copy.copy(self.bngl)
            self.bngl = "".join(self.bngl_list)

    def add_action(self, action):
        # TODO: For now I'm just assuming the "action" is just a string
        # later we'll probably make BNGAction class or something that's capable
        # of generating the string
        self.bngl = self.bngl + action + "\n"
        self.export_bngl(self.bngl_file)

    def export_bngl(self, bngl_path="model.bngl"):
        with open(bngl_path, 'w') as f:
            f.write(self.bngl)
        return bngl_path

    def find_and_set_bngl(self):
        # Let's see if we can divine the BNGL file first
        bngl_files = list(filter(lambda x: x.endswith(".bngl"), os.listdir(os.getcwd())))
        if len(bngl_files) == 0:
            print("No BNGL file given or exists in the PWD, simulator won't run")
            self.bngl_file = None
        elif len(bngl_files) > 1:
            print("More than one BNGL file exists in PWD, using {}".format(bngl_files[0]))
            self.bngl_file = bngl_files[0]
        else:
            print("BNGL file found in PWD: {}".format(bngl_files[0]))
            self.bngl_file = bngl_files[0]

    def add_params_to_bngl(self):
        self.add_action("generate_network({overwrite=>1})")
        simulate_cmd = 'simulate({'
        opts_len = len(self.run_params.keys())
        for iopt, opt in enumerate(self.run_params.keys()):
            if opt == "method":
                simulate_cmd += '{}=>"{}"'.format(opt, self.run_params[opt])
            else:
                simulate_cmd += '{}=>{}'.format(opt, self.run_params[opt])
            if iopt < opts_len-1:
                simulate_cmd += ","
        simulate_cmd += '})'
        self.add_action(simulate_cmd)
