import re, functools, subprocess, os, xmltodict, sys, shutil
import BNGUtils
from ModelStructs import Parameters, Species, MoleculeTypes, Observables, Functions,Compartments, Rules

###### CORE OBJECT AND PARSING FRONT-END ######
class BNGModel:
    '''
    The full model
    '''
    def __init__(self, bngl_model, BNGPATH=None, BNGLmode=False):
        self.active_blocks = []
        # We want blocks to be printed in the same order
        # every time
        self._action_list = ["generate_network(", "generate_hybrid_model(","simulate(", "simulate_ode(", "simulate_ssa(", "simulate_pla(", "simulate_nf(", "parameter_scan(", "bifurcate(", "readFile(", "writeFile(", "writeModel(", "writeNetwork(", "writeXML(", "writeSBML(", "writeMfile(", "writeMexfile(", "writeMDL(", "visualize(", "setConcentration(", "addConcentration(", "saveConcentration(", "resetConcentrations(", "setParameter(", "saveParameters(", "resetParameters(", "quit(", "setModelName(", "substanceUnits(", "version(", "setOption("]
        self.block_order = ["parameters", "compartments", "moltypes", 
                            "species", "observables", "functions", "rules"]
        self.BNGLmode = BNGLmode
        BNGPATH, bngexec = BNGUtils.find_BNG_path(BNGPATH)
        self.BNGPATH = BNGPATH
        self.bngexec = bngexec 
        self.model_name = ""
        self.parse_model(bngl_model)

    def __str__(self):
        '''
        write the model to str
        '''
        model_str = "begin model\n"
        for block in self.block_order:
            if block in self.active_blocks:
                model_str += str(getattr(self, block))
        model_str += "\nend model"
        return model_str

    def __repr__(self):
        return self.model_name

    def __iter__(self):
        active_ordered_blocks = [getattr(self,i) for i in self.block_order if i in self.active_blocks]
        return active_ordered_blocks.__iter__()

    def parse_model(self, model_file):
        if self.BNGLmode and model_file.endswith(".bngl"):
            # forces the old code path that tries to
            # parse bngl
            self.parse_bngl(model_file)
        else:
            # this route runs BNG2.pl on the bngl and parses
            # the XML instead
            if model_file.endswith(".bngl"):
                # TODO: Strip actions into a temp file
                # then run the gen xml 
                model_file = self.generate_xml(model_file)
                if model_file is not None:
                    self.parse_xml(model_file)
                else:
                    self.parse_bngl(model_file)
            elif model_file.endswith(".xml"):
                self.parse_xml(model_file)
            else:
                print("The extension of {} is not supported".format(model_file))
                raise NotImplemented

    def generate_xml(self, model_file):
        # TODO: Better file handling with tempfiles and 
        # not having "_stripped" in name etc
        stripped_bngl = self.strip_actions(model_file)
        rc = subprocess.run([self.bngexec, "--xml", stripped_bngl])
        if rc.returncode == 1:
            print("XML generation failed, trying the fallback parser")
            return None
        else:
            # we should now have the XML file 
            _, model_name = os.path.split(stripped_bngl)
            model_name = model_name.replace(".bngl", "")
            xml_file = model_name + ".xml"
            move_to = model_name.replace("_stripped", "") + ".xml" 
            shutil.move(xml_file, move_to)
            return move_to 

    def strip_actions(self, model_path):
        # Get model name and setup path stuff
        path, model_file = os.path.split(model_path)
        model_name = model_file.replace(".bngl","")
        stripped_file = model_name + "_stripped.bngl"
        # open model and strip actions
        with open(model_path, 'r') as mf:
            # read and strip actions
            mlines = mf.readlines()
            stripped_lines = filter(lambda x: self._not_action(x), mlines)
        # open new file and write just the model
        with open(os.path.join(path, stripped_file), 'w') as sf:
            sf.writelines(stripped_lines)
        return stripped_file

    def _not_action(self, line):
        for action in self._action_list:
            if action in line:
                return False
        return True

    def parse_xml(self, model_file):
        with open(model_file, "r") as f:
            xml_str = "".join(f.readlines())
        xml_dict = xmltodict.parse(xml_str)
        xml_model = xml_dict['sbml']['model']
        self.model_name = xml_model['@id']
        for listkey in xml_model.keys():
            if listkey == "ListOfParameters":
                param_list = xml_model[listkey]
                if param_list is not None:
                    params = param_list['Parameter']
                    self.parameters = Parameters()
                    self.parameters.parse_xml_block(params)
                    self.active_blocks.append("parameters")
            elif listkey == "ListOfObservables":
                obs_list = xml_model[listkey]
                if obs_list is not None:
                    obs = obs_list['Observable']
                    self.observables = Observables()
                    self.observables.parse_xml_block(obs)
                    self.active_blocks.append("observables")
            elif listkey == "ListOfCompartments":
                comp_list = xml_model[listkey]
                if comp_list is not None:
                    self.compartments = Compartments()
                    comps = comp_list['compartment']
                    self.compartments.parse_xml_block(comps)
                    self.active_blocks.append("compartments")
            elif listkey == "ListOfMoleculeTypes":
                mtypes_list = xml_model[listkey]
                if mtypes_list is not None:
                    mtypes = mtypes_list["MoleculeType"]
                    self.moltypes = MoleculeTypes()
                    self.moltypes.parse_xml_block(mtypes)
                    self.active_blocks.append("moltypes")
            elif listkey == "ListOfSpecies":
                species_list = xml_model[listkey]
                if species_list is not None:
                    species = species_list["Species"]
                    self.species = Species()
                    self.species.parse_xml_block(species)
                    self.active_blocks.append("species")
            elif listkey == "ListOfReactionRules":
                rrules_list = xml_model[listkey]
                if rrules_list is not None:
                    rrules = rrules_list["ReactionRule"]
                    self.rules = Rules()
                    self.rules.parse_xml_block(rrules)
                    self.active_blocks.append("rules")
            elif listkey == "ListOfFunctions":
                # TODO: Optional expression parsing?
                # TODO: Add arguments correctly
                func_list = xml_model[listkey]
                if func_list is not None:
                    self.functions = Functions()
                    funcs = func_list['Function']
                    self.functions.parse_xml_block(funcs)
                    self.active_blocks.append("functions")
        # And that's the end of parsing

    def parse_bngl(self, bngl_file):
        '''
        very basic and incomplete direct BNGL parsing
        TODO: complete a basic parsing of all blocks
        '''
        with open(bngl_file, 'r') as bngl:
            bngl_lines = bngl.readlines()

        blocks = {}
        blocks["actions"] = []
        # getting all blocks
        # TODO: handle all possible commands and options
        for iline, line in enumerate(bngl_lines):
            # TODO: Add a if statement that skips empty 
            # lines and removes it from the block line list
            if re.match(r'^begin +parameters *(#|$)', line):
                i_parameter_start = iline
            elif re.match(r'^end +parameters *(#|$)', line):
                i_parameter_end = iline
                blocks['parameters'] = bngl_lines[i_parameter_start+1:i_parameter_end]
            elif re.match(r'^begin +observables *(#|$)', line):
                i_obs_start = iline
            elif re.match(r'^end +observables *(#|$)', line):
                i_obs_end = iline
                blocks['observables'] = bngl_lines[i_obs_start+1:i_obs_end]
            elif re.match(r'^begin +compartments *(#|$)', line):
                i_compart_start = iline
            elif re.match(r'^end +compartments *(#|$)', line):
                i_compart_end = iline
                blocks['compartments'] = bngl_lines[i_compart_start+1:i_compart_end]
            elif re.match(r'^begin +molecule types *(#|$)', line):
                i_moltypes_start = iline
            elif re.match(r'^end +molecule types *(#|$)', line):
                i_moltypes_end = iline
                blocks['moltypes'] = bngl_lines[i_moltypes_start+1:i_moltypes_end]
            elif re.match(r'^begin +species *(#|$)', line):
                i_species_start = iline
            elif re.match(r'^end +species *(#|$)', line):
                i_species_end = iline
                blocks['species'] = bngl_lines[i_species_start+1:i_species_end]
            elif re.match(r'^begin +reaction +rules *(#|$)', line):
                i_rrules_start = iline
            elif re.match(r'^end +reaction +rules *(#|$)', line):
                i_rrules_end = iline
                blocks["rrules"] = bngl_lines[i_rrules_start+1:i_rrules_end]
            elif re.match(r'^begin +functions *(#|$)', line):
                i_rrules_start = iline
            elif re.match(r'^end +funcions *(#|$)', line):
                i_rrules_end = iline
                blocks["functions"] = bngl_lines[i_rrules_start+1:i_rrules_end]
            # TODO: Add other actions
            elif re.match(r'^generate_network *', line):
                blocks["actions"].append(bngl_lines[iline])
            elif re.match(r'^simulate*', line):
                blocks["actions"].append(bngl_lines[iline])
        # parse blocks
        self.parse_bngl_blocks(blocks)

    def parse_bngl_blocks(self, blocks):
        # parameters, observables, compartments, species, 
        # gen_network, simulate, rrules, functions
        for key, value in blocks.items():
            # get appropriate function to parse the block
            getattr(self, "_parse_"+key)(value)

    def _parse_parameters(self, block):
        # initialize the block object
        self.parameters = Parameters()
        # get the block to string
        self.parameters.parse_block(block)
        # active blocks list
        self.active_blocks.append("parameters")

    def _parse_species(self, block):
        # init
        self.species = Species()
        # parse
        self.species.parse_block(block)
        # add to active list
        self.active_blocks.append("species")

    def _parse_moltypes(self, block):
        # init
        self.moltypes = MoleculeTypes()
        # parse
        self.moltypes.parse_block(block)
        # add to active list
        self.active_blocks.append("moltypes")

    def _parse_observables(self, block):
        # init
        self.observables = Observables()
        # parse
        self.observables.parse_block(block)
        # add to active list
        self.active_blocks.append("observables")

    def _parse_functions(self, block):
        # init
        self.functions = Functions()
        # parse
        self.functions.parse_block(block)
        # add to active list
        self.active_blocks.append("functions")

    def _parse_rrules(self, block):
        # init
        self.rules = Rules()
        # parse
        self.rules.parse_block(block)
        # add to active list
        self.active_blocks.append("rules")

    def _parse_actions(self, block):
        # TODO: Finish this
        # print("actions")
        # print(block)
        pass

    def _parse_compartments(self, block):
        # TODO: Finish this
        compartments = list(map(lambda x: x.split(), block))
        # print("compartments")
        # print(compartments)
        pass

    def write_model(self, file_name):
        '''
        write the model to file 
        '''
        model_str = ""
        for block in self.active_blocks:
            model_str += str(getattr(self, block))
        with open(file_name, 'w') as f:
            f.write(model_str)
###### CORE OBJECT AND PARSING FRONT-END ######

if __name__ == "__main__":
    # model = BNGModel("test.bngl")
    # import IPython
    # IPython.embed()
    # with open("test.bngl", 'w') as f:
    #     f.write(str(model))
    os.chdir("validation")
    bngl_list = os.listdir(os.getcwd())
    bngl_list = filter(lambda x: x.endswith(".bngl"), bngl_list)
    for bngl in bngl_list:
        m = BNGModel(bngl)
        with open('test.bngl', 'w') as f:
            f.write(str(m))
        rc = subprocess.run([m.bngexec, 'test.bngl'])
        if rc.returncode == 1:
            print("issues with the written bngl")
            sys.exit()
    # with open("test_res.txt", "w") as f:
    #     for bngl in bngl_list:
    #         print("Working on {}".format(bngl))
    #         try:
    #             m = BNGModel(bngl)
    #         except:
    #             f.write(("Failed at {}\n".format(bngl)))
    #             # IPython.embed()
