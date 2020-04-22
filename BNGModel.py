import re, functools, subprocess, os, xmltodict, sys
import BNGUtils
from XMLPatterns import ObsPattern, MolTypePattern, RulePattern
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
        # TODO: Make sure to delete all actions, this will
        # run the model FIRST and then pull in the XML. Can 
        # lead to a ton of generated species if there is are 
        # commands in there
        rc = subprocess.run([self.bngexec, "--xml", model_file])
        if rc.returncode == 1:
            print("SBML generation failed, trying the fallback parser")
            return None
        else:
            # we should now have the SBML file 
            _, model_name = os.path.split(model_file)
            model_name = model_name.replace(".bngl", "")
            xml_file = model_name + ".xml"
            return xml_file

    def parse_xml(self, model_file):
        with open(model_file, "r") as f:
            xml_str = "".join(f.readlines())
        xml_dict = xmltodict.parse(xml_str)
        xml_model = xml_dict['sbml']['model']
        self.model_name = xml_model['@id']
        for listkey in xml_model.keys():
            if listkey == "ListOfParameters":
                param_list = xml_model[listkey]['Parameter']
                self.parameters = Parameters()
                if isinstance(param_list, list):
                    for pd in param_list:
                        self.parameters.add_item((pd['@id'],pd['@value']))
                else:
                    self.parameters.add_item((param_list['@id'], param_list['@value']))
                self.active_blocks.append("parameters")
            elif listkey == "ListOfObservables":
                obs_list = xml_model[listkey]['Observable']
                self.observables = Observables()
                # we need to turn the patterns into strings
                if isinstance(obs_list, list):
                    for od in obs_list:
                        pattern = ObsPattern(od['ListOfPatterns'])
                        self.observables.add_item((od['@type'], od['@name'], pattern))
                else: 
                    pattern = ObsPattern(obs_list)
                    self.observables.add_item((obs_list['@type'], obs_list['@name'], pattern))

                self.active_blocks.append("observables")
            elif listkey == "ListOfCompartments":
                comp_list = xml_model[listkey]
                if comp_list is not None:
                    self.compartments = Compartments()
                    comps = comp_list['compartment']
                    if isinstance(comps, list):
                        for comp in comps:
                            cname = comp['@id']
                            dim = comp['@spatialDimensions']
                            size = comp['@size']
                            if '@outside' in comp:
                                outside = comp['@outside']
                            else:
                                outside = None
                            self.compartments.add_item( (cname, dim, size, outside) )
                    else:
                        cname = comp['@id']
                        dim = comp['@spatialDimensions']
                        size = comp['@size']
                        if '@outside' in comp:
                            outside = comp['@outside']
                        else:
                            outside = None
                        self.compartments.add_item( (cname, dim, size, outside) )
                    self.active_blocks.append("compartments")
            elif listkey == "ListOfMoleculeTypes":
                mtypes_list = xml_model[listkey]["MoleculeType"]
                self.moltypes = MoleculeTypes()
                if isinstance(mtypes_list, list):
                    for md in mtypes_list:
                        pattern = MolTypePattern(md)
                        self.moltypes.add_item((pattern,))
                else:
                    pattern = MolTypePattern(md)
                    self.moltypes.add_item((pattern,))
                self.active_blocks.append("moltypes")
            elif listkey == "ListOfSpecies":
                species_list = xml_model[listkey]["Species"]
                self.species = Species()
                #TODO: Eventually regenerate patterns 
                # with bond handling from XML instead of 
                # reading the name directly to stay consistent
                for sd in species_list:
                    self.species.add_item((sd['@name'],sd['@concentration']))
                self.active_blocks.append("species")
            elif listkey == "ListOfReactionRules":
                rrules_list = xml_model[listkey]["ReactionRule"]
                self.rules = Rules()
                for rd in rrules_list:
                    rpattern = RulePattern(rd)
                    self.rules.add_item(rpattern.item_tuple)
                self.rules.consolidate_rules()
                self.active_blocks.append("rules")
            elif listkey == "ListOfFunctions":
                # TODO: Optional expression parsing?
                # TODO: Add arguments correctly
                func_list = xml_model[listkey]
                if func_list is not None:
                    self.functions = Functions()
                    funcs = func_list['Function']
                    if isinstance(funcs, list):
                         for func in funcs:
                             self.functions.add_item((func['@id'],func['Expression']))
                    else:
                         self.functions.add_item((funcs['@id'],funcs['Expression']))
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
    # model = BNGModel("validation/FceRI_ji.bngl")
    # model = BNGModel("FceRI_ji.xml")
    # model = BNGModel("egfr_net.bngl")
    model = BNGModel("egfr_net.xml")
    # model = BNGModel("compart.bngl")
    # model = BNGModel("func.bngl")
    import IPython
    IPython.embed()
    # with open("test.bngl", 'w') as f:
    #     f.write(str(model))
    # os.chdir("validation")
    # bngl_list = os.listdir(os.getcwd())
    # bngl_list = filter(lambda x: x.endswith(".bngl"), bngl_list)
    # with open("test_res.txt", "w") as f:
    #     for bngl in bngl_list:
    #         print("Working on {}".format(bngl))
    #         try:
    #             m = BNGModel(bngl)
    #         except:
    #             f.write(("Failed at {}\n".format(bngl)))
                # IPython.embed()
