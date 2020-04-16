import re
import IPython

# Objects in the model
class ModelBlock:
    def __init__(self):
        self._item_list = []

    def add_item(self, item_tpl):
        # TODO: try adding evaluation of the parameter here
        # for the future, in case we want people to be able
        # to adjust the math
        # TODO: Error handling, some names will definitely break this
        name, value = item_tpl
        try:
            setattr(self, name, value)
        except:
            pass
        self._item_list.append((name, value))

    def add_items(self, item_list):
        for item in item_list:
            self.add_item(item)

    def print_all(self):
        print(self._item_list)

    def write_block(self):
        block_lines = ["begin {}".format(self.name)]
        for item in self._item_list:
            block_lines.append("  " + " ".join(item))
        block_lines.append("end {}".format(self.name))
        return block_lines


# TODO: Add a LOT of error handling
# TODO: Each object requires a write function 
# appropriate to the object
class Parameters(ModelBlock):
    '''
    Class containing parameters
    '''
    def __init__(self):
        super().__init__()
        self.name = "parameters"

class Species(ModelBlock):
    '''
    Class containing species
    '''
    def __init__(self):
        super().__init__()
        self.name = "species"

class MoleculeTypes(ModelBlock):
    '''
    Class containing molecule types 
    '''
    def __init__(self):
        super().__init__()
        self.name = "molecule types"

    def add_item(self, name):
        self._item_list.append(tuple(name))

class Observables(ModelBlock):
    '''
    Class for observables
    '''
    def __init__(self):
        super().__init__()
        self.name = "species"

    def add_item(self, item_tpl): 
        otype, name, pattern = item_tpl
        self._item_list.append((otype, name, pattern))

class Functions(ModelBlock):
    '''
    Class for functions
    '''
    def __init__(self):
        super().__init__()
        self.name = "functions"

class Rules(ModelBlock):
    def __init__(self):
        super().__init__()
        self.name = "reaction rules"

    def add_item(self, item_tpl):
        rule_txt = item_tpl
        self._item_list.append(rule_txt)

# Now onto the actual model and parsing
class BNGModel:
    '''
    The full model
    '''
    def __init__(self, bngl_file):
        self.changed = False
        self.active_blocks = []
        self.parse_bngl(bngl_file)
        self.write_model()

    def parse_bngl(self, bngl_file):
        with open(bngl_file, 'r') as bngl:
            bngl_lines = bngl.readlines()

        blocks = {}
        blocks["actions"] = []
        # getting all blocks
        # TODO: handle all possible commands and options
        # options currently are:
        # parameters, observables, compartments, species, 
        # gen_network, simulate, rrules, molecule types
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
        self.parse_blocks(blocks)

    def parse_blocks(self, blocks):
        # parameters, observables, compartments, species, 
        # gen_network, simulate, rrules, functions
        for key, value in blocks.items():
            # get appropriate function to parse the block
            getattr(self, "_parse_"+key)(value)

    def strip_comment(self, line):
        return line[0:line.find("#")]

    def _parse_parameters(self, block):
        # strip comments
        params = list(map(self.strip_comment, block))
        # split properly
        params = list(map(lambda x: x.split(), params))
        self.parameters = Parameters()
        self.parameters.add_items(params)
        print("parameters")
        self.parameters.print_all()
        self.active_blocks.append("parameters")

    def _parse_species(self, block):
        # strip comments 
        species = list(map(self.strip_comment, block))
        # basically the same as parameters 
        species = list(map(lambda x: x.split(), species))
        self.species = Species()
        self.species.add_items(species)
        print("species")
        self.species.print_all()

    def _parse_moltypes(self, block):
        # strip comments 
        moltypes = list(map(self.strip_comment, block))
        # remove white spaces 
        moltypes = list(map(lambda x: x.split(), moltypes))
        self.moltypes = MoleculeTypes()
        self.moltypes.add_items(moltypes)
        print("molecule types")
        self.moltypes.print_all()

    def _parse_observables(self, block):
        # strip comments 
        obs = list(map(self.strip_comment, block))
        # remove white spaces and split
        obs = list(map(lambda x: x.split(), obs))
        self.observables = Observables()
        self.observables.add_items(obs)
        print("observables")
        self.observables.print_all()

    def _parse_functions(self, block):
        # strip comments
        functions = list(map(self.strip_comment, block))
        # split by = sign
        functions = list(map(lambda x: x.split("="), functions))
        self.functions = Functions()
        self.functions.add_items(functions)
        self.functions.print_all()

    def _parse_rrules(self, block):
        # strip comments
        rules = list(map(self.strip_comment, block))
        # split by = sign
        rules = list(map(lambda x: " ".join(x.split()), rules))
        self.rules = Rules()
        self.rules.add_items(rules)
        self.rules.print_all()

    def _parse_actions(self, block):
        print("actions")
        print(block)

    def _parse_compartments(self, block):
        compartments = list(map(lambda x: x.split(), block))
        print("compartments")
        print(compartments)

    def write_model(self):
        '''
        write the model to str
        '''
        model_str_list = []
        for block in self.active_blocks:
            model_str_list += getattr(self, block).write_block()
        print("\n".join(model_str_list))
