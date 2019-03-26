import os

class BNGParser:
    '''
    This is a very simple BNGL reader that's able to read,y write the file, 
    strip actions that are in the file originally in there and able to append
    lines to the file. 
    '''
    def __init__(self, bngl):
        if bngl is None:
            self.find_and_set_bngl()
        else:
            self.bngl_file = bngl
        if self.bngl_file is not None:
            self.init_bngl(self.bngl_file)

    def init_bngl(self, bngl_file):
        self.bngl_list = self.read_to_list(bngl_file)
        self.bngl = "".join(self.bngl_list)

    def read_to_list(self, bngl_file):
        with open(bngl_file, 'r') as f:
            bngl_list = f.readlines()
        return bngl_list

    def clean_actions(self):
        no_action_bngl = []
        for line in self.bngl_list:
            if "end model" in line: 
                no_action_bngl.append(line)
                break
            no_action_bngl.append(line)
        self.bngl_list = no_action_bngl
        self.bngl = "".join(self.bngl_list)

    def add_action(self, action):
        # TODO: For now I'm just assuming the "action" is just a string
        # later we'll probably make BNGAction class or something that's capable
        # of generating the string
        self.bngl = self.bngl + action + "\n"

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
