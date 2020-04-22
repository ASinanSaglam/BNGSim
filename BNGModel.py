import re, functools, subprocess, os, xmltodict, sys
import BNGUtils
import IPython

###### BONDS #####
class Bonds:
    def __init__(self, bonds_xml=None):
        self.bonds_dict = {}
        if bonds_xml is not None:
            self.resolve_xml(bonds_xml)

    def set_xml(self, bonds_xml):
        self.resolve_xml(bonds_xml)

    def get_bond_id(self, comp):
        # Get the ID of the bond from an XML id something 
        # belongs to, e.g. O1_P1_M1_C2 
        num_bonds = comp["@numberOfBonds"]
        comp_id = comp["@id"]
        try: 
            num_bond = int(num_bonds)
        except: 
            # This means we have something like +/?
            return num_bonds
        # use the comp_id to find the bond index from 
        # self.bonds_dict 
        comp_key = self.get_tpl_from_id(comp_id)
        bond_id = self.bonds_dict[comp_key]
        return bond_id
        
    def get_tpl_from_id(self, id_str):
        # ID str is looking like O1_P1_M2_C3
        # we are going to assume a 4-tuple per key
        id_list = id_str.split("_")
        # id_tpl = tuple([self.get_tpl(x) for x in id_list])
        id_tpl = tuple(id_list)
        return id_tpl

    def tpls_from_bond(self, bond):
        s1 = bond["@site1"] 
        s2 = bond["@site2"]
        id_list_1 = s1.split("_")
        #s1_tpl = tuple([int(x[1:]) for x in id_list_1])
        s1_tpl = tuple(id_list_1)
        id_list_2 = s2.split("_")
        #s2_tpl = tuple([int(x[1:]) for x in id_list_2])
        s2_tpl = tuple(id_list_2)
        return (s1_tpl, s2_tpl) 

    def resolve_xml(self, bonds_xml):
        # OK, we want to unfuck this code. 
        # self.bonds_dict is a dictionary you can key
        # with the tuple taken from the ID and then 
        # get a bond ID cleanly
        if isinstance(bonds_xml, list):
            for ibond, bond in enumerate(bonds_xml): 
                bond_partner_1, bond_partner_2 = self.tpls_from_bond(bond)
                self.bonds_dict[bond_partner_1] = ibond+1
                self.bonds_dict[bond_partner_2] = ibond+1
        else:
            bond_partner_1, bond_partner_2 = self.tpls_from_bond(bonds_xml)
            self.bonds_dict[bond_partner_1] = 1
            self.bonds_dict[bond_partner_2] = 1
###### BONDS #####

###### PATTERNS ###### 
class Pattern:
    def __init__(self, pattern_xml):
        self.pattern_xml = pattern_xml
        self.bonds = Bonds()
        self.string = self.resolve_xml(self.pattern_xml)

    def ___repr__(self):
        return self.string

    def __str__(self):
        return self.string

    def mol_to_str(self, mol_xml):
        # TODO: Document, especially the bond mechanism
        if not isinstance(mol_xml, list):
            mol_str = mol_xml["@name"] + "("
            if "ListOfComponents" in mol_xml:
                # Single molecule can't have bonds
                comp_str = self.comp_to_str(mol_xml["ListOfComponents"]["Component"])
                mol_str += comp_str
            mol_str += ")"
        else:
            # this means we have multiple molecules bonding
            mol_str = ""
            bonds_to_pass = {}
            for imol, mol in enumerate(mol_xml):
                if imol > 0:
                    # complexing
                    mol_str += "."
                mol_str += mol["@name"] + "("
                if "ListOfComponents" in mol:
                    comp_str = self.comp_to_str(mol['ListOfComponents']['Component'])
                    mol_str += comp_str
                mol_str += ")"
        return mol_str

    def comp_to_str(self, comp_xml):
        # bonds = compartment id, bond id 
        # comp xml can be a list or a dict
        if isinstance(comp_xml, list):
            # we have multiple and this is a list
            comp_str = ""
            for icomp, comp in enumerate(comp_xml):
                comp_id = int(comp["@id"].split("_")[-1].replace("C",""))
                if icomp > 0:
                    comp_str += ","
                comp_str += comp['@name']
                if "@state" in comp:
                    comp_str += "~{}".format(comp['@state'])
                if comp["@numberOfBonds"] != '0':
                    bond_id = self.bonds.get_bond_id(comp)
                    comp_str += "!{}".format(bond_id)
        else:
            # single comp, this is a dict
            comp_str = comp_xml['@name']
            comp_id = int(comp_xml["@id"].split("_")[-1].replace("C",""))
            if "@state" in comp_xml:
                comp_str += "~{}".format(comp_xml['@state'])
            if comp_xml['@numberOfBonds'] != '0':
                bond_id = self.bonds.get_bond_id(comp_xml)
                comp_str += "!{}".format(bond_id)
        return comp_str

class ObsPattern(Pattern):
    def __init__(self, pattern_xml):
        super().__init__(pattern_xml)

    def resolve_xml(self, obs_xml):
        patterns = obs_xml['Pattern']
        if isinstance(patterns, list):
            # we have multiple stuff so this becomes a list
            obs_str = ""
            for ipattern, pattern in enumerate(patterns): 
                # 
                if "ListOfBonds" in pattern:
                    self.bonds.set_xml(pattern["ListOfBonds"]["Bond"])
                if ipattern > 0:
                    obs_str += ","
                mol = pattern['ListOfMolecules']['Molecule']
                obs_res = self.mol_to_str(mol)
                obs_str += obs_res
        else:
            bonds_list = []
            if "ListOfBonds" in patterns:
                self.bonds.set_xml(patterns["ListOfBonds"]["Bond"])
            mol = patterns['ListOfMolecules']["Molecule"]
            obs_str = self.mol_to_str(mol)
        return obs_str

class MolTypePattern(Pattern):
    def __init__(self, pattern_xml):
        super().__init__(pattern_xml)

    def resolve_xml(self, xml_molt):
        molt_str = xml_molt['@id'] + "("
        comp_dict = xml_molt['ListOfComponentTypes']['ComponentType']
        if '@id' in comp_dict:
            molt_str += comp_dict['@id']
        else:
            # multiple components
            for icomp, comp in enumerate(comp_dict):
                if icomp > 0:
                    molt_str += ","
                molt_str += comp['@id']
                if "ListOfAllowedStates" in comp:
                    # we have states
                    al_states = comp['ListOfAllowedStates']['AllowedState']
                    for istate, state in enumerate(al_states):
                        molt_str += "~{}".format(state['@id'])
        molt_str += ")"
        return molt_str

class RulePattern(Pattern):
    # TODO: Add bond handling
    def __init__(self, pattern_xml):
        super().__init__(pattern_xml)

    def resolve_xml(self, pattern_xml):
        '''
        in this particular case also sets the self.item_tuple
        '''
        # 
        rule_name = pattern_xml['@name']
        lhs = self.resolve_rxn_side(pattern_xml['ListOfReactantPatterns'])
        rhs = self.resolve_rxn_side(pattern_xml['ListOfProductPatterns'])
        rate_law = self.resolve_ratelaw(pattern_xml['RateLaw'])
        # 
        # We need to set self.item_tuple to 
        # (rule_name, LHS, RHS, rate_law)
        self.item_tuple = (rule_name, lhs, rhs, rate_law)
        return "{}: {} -> {} {}".format(rule_name, lhs, rhs, rate_law)

    def resolve_ratelaw(self, rate_xml):
        rate_type = rate_xml['@type']
        rate_cts_xml = rate_xml['ListOfRateConstants']
        # TODO: what happens if this is not just a simple 
        # rate constant? Probably this becomes a list 
        # and we'll fail right here.
        rate_cts = rate_cts_xml['RateConstant']['@value']
        return rate_cts

    def resolve_rxn_side(self, side_xml):
        # this is either reactant or product
        if 'ReactantPattern' in side_xml:
            # this is a lhs/reactant side
            side_list = side_xml['ReactantPattern']
            if isinstance(side_list, list):
                # this is a list of reactants
                react_str = ""
                for ireact, react in enumerate(side_list):
                    if "ListOfBonds" in react:
                        self.bonds.set_xml(react["ListOfBonds"]['Bond'])
                    if ireact > 0:
                        react_str += " + "
                    # bonds should go here
                    react_res = self.mol_to_str(react['ListOfMolecules']['Molecule'])
                    react_str += react_res
            else: 
                if "ListOfBonds" in side_list:
                    self.bonds.set_xml(side_list["ListOfBonds"]['Bond'])
                # TODO: a single item? 
                react_res = self.mol_to_str(side_list['ListOfMolecules']['Molecule'])
                react_str = react_res
            return react_str
        elif "ProductPattern" in side_xml:
            side_list = side_xml['ProductPattern']
            if isinstance(side_list, list):
                # this is a list of reactants
                prod_str = ""
                for iprod, prod in enumerate(side_list):
                    if "ListOfBonds" in prod:
                        self.bonds.set_xml(prod["ListOfBonds"]['Bond'])
                    if iprod > 0:
                        prod_str += " + "
                    prod_res = self.mol_to_str(prod['ListOfMolecules']['Molecule'])
                    prod_str += prod_res
            else: 
                # TODO: a single item? 
                if "ListOfBonds" in side_list:
                    self.bonds.set_xml(side_list["ListOfBonds"]['Bond'])
                prod_res = self.mol_to_str(side_list['ListOfMolecules']['Molecule'])
                prod_str = prod_res
            return prod_str
        else: 
            print("do not recognize XML: {}".format(side_xml))
###### PATTERNS ###### 

###### MODEL STRUCTURES ###### 
# Objects in the model
class ModelBlock:
    def __init__(self):
        self._item_dict = {}

    def __len__(self):
        return len(self._item_dict.keys())

    def __repr__(self):
        # overwrites what the class representation
        # shows the items in the model block in 
        # say ipython
        return str(self._item_dict)

    def add_item(self, item_tpl):
        # TODO: try adding evaluation of the parameter here
        # for the future, in case we want people to be able
        # to adjust the math
        # TODO: Error handling, some names will definitely break this
        name, value = item_tpl
        self._item_dict[name] = value
        try:
            setattr(self, name, value)
        except:
            print("can't set {} to {}".format(name, value))
            pass

    def add_items(self, item_list):
        for item in item_list:
            self.add_item(item)

    def print(self):
        print(self)

    def strip_comment(self, line):
        return line[0:line.find("#")]

    def parse_block(self, block):
        raise NotImplemented


# TODO: Add a LOT of error handling
class Parameters(ModelBlock):
    '''
    Class containing parameters
    '''
    def __init__(self):
        super().__init__()
        self.name = "parameters"

    def __setattr__(self, name, value):
        if hasattr(self, "_item_dict"):
            if name in self._item_dict.keys():
                self._item_dict[name] = value
        self.__dict__[name] = value

    def __str__(self):
        # overwrites what the method returns when 
        # it's converted to string
        block_lines = ["\nbegin {}".format(self.name)]
        for item in self._item_dict.keys():
            block_lines.append("  " + "{} {}".format(item, self._item_dict[item]))
        block_lines.append("end {}\n".format(self.name))
        return "\n".join(block_lines)

    def parse_block(self, block):
        # strip comments
        params = list(map(self.strip_comment, block))
        # split properly
        params = list(map(lambda x: x.split(), params))
        # ensure we have at least 2 elements
        params = list(filter(lambda x: len(x)>1, params))
        # list should be param_name = expression
        params = list(map(lambda x: [x[0], str(functools.reduce(lambda x,y: x+y, x[1:]))], params))
        # now generate the param object
        self.add_items(params)

class Species(ModelBlock):
    '''
    Class containing species
    '''
    def __init__(self):
        super().__init__()
        self.name = "species"

    def __str__(self):
        # overwrites what the method returns when 
        # it's converted to string
        block_lines = ["\nbegin {}".format(self.name)]
        for item in self._item_dict.keys():
            block_lines.append("  " + "{} {}".format(item,self._item_dict[item]))
        block_lines.append("end {}\n".format(self.name))
        return "\n".join(block_lines)

    def parse_block(self, block):
        # strip comments 
        species = list(map(self.strip_comment, block))
        # basically the same as parameters 
        species = list(map(lambda x: x.split(), species))
        self.add_items(species)

    def __getitem__(self, key):
        return self._item_dict[key]

    def __setitem__(self, key, value):
        self._item_dict[key] = value

class MoleculeTypes(ModelBlock):
    '''
    Class containing molecule types 
    '''
    def __init__(self):
        super().__init__()
        self.name = "molecule types"

    def add_item(self, item_tpl):
        name, = item_tpl
        self._item_dict[name] = ""

    def __str__(self):
        # overwrites what the method returns when 
        # it's converted to string
        block_lines = ["\nbegin {}".format(self.name)]
        for item in self._item_dict.keys():
            block_lines.append("  " + "{}".format(item))
        block_lines.append("end {}\n".format(self.name))
        return "\n".join(block_lines)

    def parse_block(self, block):
        # strip comments 
        moltypes = list(map(self.strip_comment, block))
        # remove white spaces 
        moltypes = list(map(lambda x: x.split(), moltypes))
        self.add_items(moltypes)

class Observables(ModelBlock):
    '''
    Class for observables
    '''
    def __init__(self):
        super().__init__()
        self.name = "observables"

    def __setattr__(self, name, value):
        if hasattr(self, "_item_dict"):
            if name in self._item_dict.keys():
                self._item_dict[name][1] = value
        self.__dict__[name] = value

    def add_item(self, item_tpl): 
        otype, name, pattern = item_tpl
        self._item_dict[name] = [otype, pattern]
        try:
            setattr(self, name, pattern)
        except:
            print("can't set {} to {}".format(name, pattern))
            pass

    def __str__(self):
        # overwrites what the method returns when 
        # it's converted to string
        block_lines = ["\nbegin {}".format(self.name)]
        for item in self._item_dict.keys():
            block_lines.append("  " + 
                    "{} {} {}".format(self._item_dict[item][0],
                                      item,
                                      self._item_dict[item][1]))
        block_lines.append("end {}\n".format(self.name))
        return "\n".join(block_lines)

    def parse_block(self, block):
        # strip comments 
        obs = list(map(self.strip_comment, block))
        # remove white spaces and split
        obs = list(map(lambda x: x.split(), obs))
        self.add_items(obs)


class Functions(ModelBlock):
    '''
    Class for functions
    '''
    def __init__(self):
        super().__init__()
        self.name = "functions"

    # TODO: Fix this such that we can re-write functions
    def __str__(self):
        # overwrites what the method returns when 
        # it's converted to string
        block_lines = ["\nbegin {}".format(self.name)]
        print(self._item_dict)
        for item in self._item_dict.keys():
            block_lines.append("  " + 
                    "{} {} {}".format(self._item_dict[item][0],
                                      item,
                                      self._item_dict[item][1]))
        block_lines.append("end {}\n".format(self.name))
        return "\n".join(block_lines)

    def parse_block(self, block):
        # strip comments
        functions = list(map(self.strip_comment, block))
        # split by = sign
        functions = list(map(lambda x: x.split("="), functions))
        self.functions.add_items(functions)

class Rules(ModelBlock):
    def __init__(self):
        super().__init__()
        self.name = "reaction rules"

    def add_item(self, item_tpl):
        # TODO: handle this entirely differently and 
        # properly parse rules
        '''
        A reaction rule, requires a 5-tuple
        (rule_name, LHS, RHS, rate_law)
        '''
        rule_name, lhs, rhs, rate_law = item_tpl
        self._item_dict[rule_name] = (lhs, rhs, rate_law)

    def __str__(self):
        # TODO: printing also needs a lot of adjusting
        block_lines = ["\nbegin {}".format(self.name)]
        for item in self._item_dict.keys():
            rule_tpl = self._item_dict[item]
            rule_str = ""
            if item != "":
                rule_str += "  {}: ".format(item)
            rule_str += rule_tpl[0]
            rule_str += " -> "
            rule_str += rule_tpl[1] 
            rule_str += " {}".format(rule_tpl[2])
            block_lines.append(rule_str)
        block_lines.append("end {}\n".format(self.name))
        return "\n".join(block_lines)

    def parse_block(self, block):
        # TODO: Update this so it parses the rules to at least
        # some extent. Don't need anything complicated, just need 
        # a barebones implementation 

        # strip comments
        rules = list(map(self.strip_comment, block))
        # split 
        rules = list(map(lambda x: " ".join(x.split(" ")), rules))
        # FIXME: This needs a 4-tuple per item, see above
        self.add_items(rules)
###### MODEL STRUCTURES ###### 

###### CORE OBJECT AND PARSING FRONT-END ######
# Now onto the actual model and parsing
class BNGModel:
    '''
    The full model
    '''
    def __init__(self, bngl_model, BNGPATH=None, BNGLmode=False):
        self.active_blocks = []
        self.BNGLmode = BNGLmode
        BNGPATH, bngexec = BNGUtils.find_BNG_path(BNGPATH)
        self.BNGPATH = BNGPATH
        self.bngexec = bngexec 
        self.parse_model(bngl_model)

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
        # TODO: implement the XML parser
        print("Parsing the XML")
        with open(model_file, "r") as f:
            xml_str = "".join(f.readlines())
        xml_dict = xmltodict.parse(xml_str)
        xml_model = xml_dict['sbml']['model']
        for listkey in xml_model.keys():
            if listkey == "ListOfParameters":
                param_list = xml_model[listkey]['Parameter']
                self.parameters = Parameters()
                for pd in param_list:
                    self.parameters.add_item((pd['@id'],pd['@value']))
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
                    # TODO: implement compartment parsing
                    pass
                    self.active_blocks.append("compartments")
            elif listkey == "ListOfMoleculeTypes":
                mtypes_list = xml_model[listkey]["MoleculeType"]
                self.moltypes = MoleculeTypes()
                for md in mtypes_list:
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
                # TODO: We need to turn these into strings
                # FIXME: Temporarily passing this to 
                # fix observables
                for rd in rrules_list:
                    # we will need a 4-tuple
                    rpattern = RulePattern(rd)
                    self.rules.add_item(rpattern.item_tuple)
                self.active_blocks.append("rules")
            elif listkey == "ListOfFunctions":
                #TODO: Implement functions
                pass
        # Tons more work to do with implementing XML parsing
        # IPython.embed()

    def __str__(self):
        '''
        write the model to str
        '''
        model_str = "begin model\n"
        for block in self.active_blocks:
            model_str += str(getattr(self, block))
        model_str += "\nend model"
        return model_str

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
        # print(model_str)
        with open(file_name, 'w') as f:
            f.write(model_str)
###### CORE OBJECT AND PARSING FRONT-END ######

if __name__ == "__main__":
    # model = BNGModel("validation/FceRI_ji.bngl")
    # model = BNGModel("FceRI_ji.xml")
    # model = BNGModel("egfr_net.bngl")
    model = BNGModel("egfr_net.xml")
    IPython.embed()
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
