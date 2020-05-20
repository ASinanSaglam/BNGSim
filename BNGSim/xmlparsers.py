from BNGSim.bonds import Bonds
from BNGSim.molecule import Molecule

###### PATTERNS ###### 
class XMLObj:
    def __init__(self, pattern_xml):
        self.molecules = []
        self.pattern_xml = pattern_xml
        self.bonds = Bonds()
        self.string = self.resolve_xml(self.pattern_xml)

    def __repr__(self):
        return self.gen_string()

    def __str__(self):
        return self.gen_string()

    def gen_string(self):
        # TODO: Eventually slap in a mechanism where
        # this only runs when the underlying structure 
        # changes 
        pat_str = ""
        for imol, mol in enumerate(self.molecules):
            if imol > 0:
                pat_str += "."
            pat_str += str(mol)
        return pat_str

    def mol_to_str(self, mol_xml):
        # TODO: Document, especially the bond mechanism
        mol_list = []
        if not isinstance(mol_xml, list):
            # we are going to store molecules, components
            # and compartments in a separate dictionary 
            # for use later
            name = mol_xml['@name'] 
            # molecule dictionary
            mol_dict = {"name": name, "components": [], "compartment": None, "outer_comp": None} 
            # start building the string
            mol_str = mol_xml["@name"] 
            if "ListOfComponents" in mol_xml:
                mol_str += "("
                # Single molecule can't have bonds
                comp_str = self.comp_to_str(mol_xml["ListOfComponents"]["Component"], mol_dict=mol_dict)
                mol_str += comp_str
                mol_str += ")"
            if '@compartment' in mol_xml:
                mol_dict['compartment'] = mol_xml['@compartment']
                mol_str += "@{}".format(mol_xml['@compartment'])
            mol_list.append(mol_dict)
        else:
            # this means we have multiple molecules bonding
            mol_str = ""
            for imol, mol in enumerate(mol_xml):
                name = mol['@name']
                # we can have multiple copies with different 
                # component states bound together. we need to
                # somehow account for it
                mol_dict = {"name": name, "components": [], "compartment": None, "outer_comp": None}
                if imol > 0:
                    # complexing
                    mol_str += "."
                mol_str += mol["@name"] 
                if "ListOfComponents" in mol:
                    mol_str += "("
                    comp_str = self.comp_to_str(mol['ListOfComponents']['Component'], mol_dict=mol_dict)
                    mol_str += comp_str
                    mol_str += ")"
                if '@compartment' in mol:
                    mol_dict['compartment'] = mol['@compartment']
                    mol_str += "@{}".format(mol['@compartment'])
                mol_list.append(mol_dict)
        mol_obj = Molecule(mol_list)
        return mol_str, mol_obj

    def comp_to_str(self, comp_xml, mol_dict=None):
        # bonds = compartment id, bond id 
        # comp xml can be a list or a dict
        if isinstance(comp_xml, list):
            # we have multiple and this is a list
            comp_str = ""
            for icomp, comp in enumerate(comp_xml):
                comp_dict = {}
                if icomp > 0:
                    comp_str += ","
                comp_str += comp['@name']
                comp_dict['name'] = comp['@name']
                if "@label" in comp:
                    comp_dict['label'] = comp['@label']
                    comp_str += "%{}".format(comp['@label'])
                if "@state" in comp:
                    comp_dict['state'] = comp['@state']
                    comp_str += "~{}".format(comp['@state'])
                if comp["@numberOfBonds"] != '0':
                    comp_dict['bonds'] = []
                    bond_id = self.bonds.get_bond_id(comp)
                    for bi in bond_id:
                        comp_dict['bonds'].append(bi)
                        comp_str += "!{}".format(bi)
                if mol_dict is not None:
                    mol_dict['components'].append(comp_dict)
        else:
            # single comp, this is a dict
            comp_dict = {}
            comp_str = comp_xml['@name']
            comp_dict['name'] = comp_xml['@name']
            if "@label" in comp_xml:
                comp_dict['label'] = comp_xml['@label']
                comp_str += "%{}".format(comp_xml['@label'])
            if "@state" in comp_xml:
                comp_dict['state'] = comp_xml['@state']
                comp_str += "~{}".format(comp_xml['@state'])
            if comp_xml['@numberOfBonds'] != '0':
                comp_dict['bonds'] = []
                bond_id = self.bonds.get_bond_id(comp_xml)
                for bi in bond_id:
                    comp_dict['bonds'].append(bi)
                    comp_str += "!{}".format(bi)
            if mol_dict is not None:
                mol_dict['components'].append(comp_dict)
        return comp_str

class ObsXML(XMLObj):
    def __init__(self, pattern_xml):
        super().__init__(pattern_xml)

    def gen_string(self):
        obs_str = ""
        for iobs, obs in enumerate(self.molecules):
            if iobs > 0:
                obs_str += ","
            obs_str += str(obs)
        return obs_str

    def resolve_xml(self, obs_xml):
        patterns = obs_xml['Pattern']
        if isinstance(patterns, list):
            # we have multiple stuff so this becomes a list
            obs_str = ""
            obs_list = []
            for ipattern, pattern in enumerate(patterns): 
                # 
                if "ListOfBonds" in pattern:
                    self.bonds.set_xml(pattern["ListOfBonds"]["Bond"])
                if ipattern > 0:
                    obs_str += ","
                mol = pattern['ListOfMolecules']['Molecule']
                obs_res, mol_obj = self.mol_to_str(mol)
                if '@compartment' in pattern:
                    obs_str += "@{}:".format(pattern['@compartment'])
                    mol_obj.set_outer_compt(pattern['@compartment'])
                self.molecules.append(mol_obj)
                obs_str += obs_res
        else:
            if "ListOfBonds" in patterns:
                self.bonds.set_xml(patterns["ListOfBonds"]["Bond"])
            mol = patterns['ListOfMolecules']["Molecule"]
            obs_str = ""
            obs_res, mol_obj = self.mol_to_str(mol)
            if '@compartment' in patterns:
                obs_str += "@{}:".format(patterns['@compartment'])
                mol_obj.set_outer_compt(patterns['@compartment'])
            obs_str += obs_res
            self.molecules.append(mol_obj)
        return obs_str

class SpeciesXML(XMLObj):
    def __init__(self, pattern_xml):
        super().__init__(pattern_xml)

    def resolve_xml(self, spec_xml):
        pattern = spec_xml['ListOfMolecules']['Molecule']
        if "@compartment" in spec_xml:
            outer_comp = spec_xml["@compartment"]
        else:
            outer_comp = None
        # bonds stored in spec_xml
        if "ListOfBonds" in spec_xml:
            self.bonds.set_xml(spec_xml["ListOfBonds"]["Bond"])
        # list of a singular species?
        mol_list = []
        if isinstance(pattern, list):
            spec_str = ""
            for ipat, pat in enumerate(pattern): 
                if ipat > 0:
                    spec_str += "."
                spec_res, ml = self.mol_to_str(pat)
                if '@compartment' in pat:
                    spec_str += "@{}:".format(pat['@compartment'])
                if ipat == 0:
                    ml.set_outer_compt(outer_comp)
                mol_list.append(ml)
                spec_str += spec_res
            self.molecules += mol_list
        else:
            spec_str, mol_obj = self.mol_to_str(pattern)
            mol_obj.set_outer_compt(outer_comp)
            self.molecules.append(mol_obj)
        return spec_str

class MolTypeXML(XMLObj):
    def __init__(self, pattern_xml):
        super().__init__(pattern_xml)

    def gen_string(self):
        return str(self.molecules[0])

    def resolve_xml(self, molt_xml):
        molt_str = molt_xml['@id'] 
        mol_list = []
        mol_dict = {"name": molt_str, "components": [], "compartment": None, "outer_comp": None}
        if 'ListOfComponentTypes' in molt_xml:
            molt_str += "("
            comp_dict = molt_xml['ListOfComponentTypes']['ComponentType']
            if '@id' in comp_dict:
                molt_str += comp_dict['@id']
                cd = {"name": comp_dict["@id"]}
                if "ListOfAllowedStates" in comp_dict:
                    # we have states
                    al_states = comp_dict['ListOfAllowedStates']['AllowedState']
                    cd["states"] = []
                    if isinstance(al_states, list):
                        for istate, state in enumerate(al_states):
                            molt_str += "~{}".format(state['@id'])
                            cd["states"].append(state["@id"])
                    else:
                        molt_str += "~{}".format(al_states['@id'])
                        cd["states"].append(state["@id"])
                mol_dict["components"].append(cd)
            else:
                # multiple components
                for icomp, comp in enumerate(comp_dict):
                    if icomp > 0:
                        molt_str += ","
                    molt_str += comp['@id']
                    cd = {"name": comp['@id']}
                    if "ListOfAllowedStates" in comp:
                        # we have states
                        al_states = comp['ListOfAllowedStates']['AllowedState']
                        cd['states'] = []
                        if isinstance(al_states, list):
                            for istate, state in enumerate(al_states):
                                molt_str += "~{}".format(state['@id'])
                                cd['states'].append(state['@id'])
                        else:
                            molt_str += "~{}".format(al_states['@id'])
                            cd['states'].append(state['@id'])
                    mol_dict['components'].append(cd)
            molt_str += ")"
        mol_list.append(mol_dict)
        self.molecules.append(Molecule(mol_list))
        return molt_str

class RuleXML(XMLObj):
    def __init__(self, pattern_xml):
        self.bidirectional = False
        super().__init__(pattern_xml)

    def gen_string(self):
        if self.bidirectional:
            return "{}: {} <-> {} {},{}".format(self.name, self.side_string(self.lhs_list), self.side_string(self.rhs_list), self.rate_law[0], self.rate_law[1])
        else:
            return "{}: {} -> {} {}".format(self.name, self.side_string(self.lhs_list), self.side_string(self.rhs_list), self.rate_law[0])

    def side_string(self, side_list):
        pat_str = ""
        for ipat, pat in enumerate(side_list):
            if ipat > 0:
                pat_str += " + "
            pat_str += str(pat)
        return pat_str

    def set_rate_law(self, rate_law):
        if len(rate_law) == 1:
            self.rate_law = [rate_law[0]]
        elif len(rate_law) == 2: 
            self.rate_law = [rate_law[0], rate_law[1]]
            self.bidirectional = True
        else:
            print("1 or 2 rate constants allowed")
    
    def resolve_xml(self, pattern_xml):
        '''
        '''
        # 
        rule_name = pattern_xml['@name']
        self.name = rule_name
        self.lhs, self.lhs_list = self.resolve_rxn_side(pattern_xml['ListOfReactantPatterns'])
        self.rhs, self.rhs_list = self.resolve_rxn_side(pattern_xml['ListOfProductPatterns'])
        if 'RateLaw' not in pattern_xml:
            print("Rule seems to be missing a rate law, please make sure that XML exporter of BNGL supports whatever you are doing!")
        self.rate_law = [self.resolve_ratelaw(pattern_xml['RateLaw'])]

    def resolve_ratelaw(self, rate_xml):
        rate_type = rate_xml['@type']
        if rate_type == 'Ele':
            rate_cts_xml = rate_xml['ListOfRateConstants']
            rate_cts = rate_cts_xml['RateConstant']['@value']
        elif rate_type == 'Function':
            rate_cts = rate_xml['@name']
        elif rate_type == 'MM' or rate_type == "Sat":
            # A function type 
            rate_cts = rate_type + "("
            args = rate_xml['ListOfRateConstants']["RateConstant"]
            if isinstance(args, list):
                for iarg, arg in enumerate(args):
                    if iarg > 0:
                        rate_cts += ","
                    rate_cts += arg["@value"]
            else:
                rate_cts += args["@value"]
            rate_cts += ")"
        else:
            print("don't recognize rate law type")
        return rate_cts

    def resolve_rxn_side(self, side_xml):
        # this is either reactant or product
        if side_xml is None:
            return "0", [[None]]
        elif 'ReactantPattern' in side_xml:
            # this is a lhs/reactant side
            sl = []
            side_list = side_xml['ReactantPattern']
            if '@compartment' in side_list:
                react_str = "@{}:".format(side_list['@compartment'])
                outer_comp = side_list['@compartment']
            else:
                react_str = ""
                outer_comp = None
            if isinstance(side_list, list):
                # side list to save
                # this is a list of reactants
                for ireact, react in enumerate(side_list):
                    if "ListOfBonds" in react:
                        self.bonds.set_xml(react["ListOfBonds"]['Bond'])
                    if ireact > 0:
                        react_str += " + "
                    # bonds should go here
                    react_res, mol_obj = self.mol_to_str(react['ListOfMolecules']['Molecule'])
                    mol_obj.set_outer_compt(outer_comp)
                    sl.append(mol_obj)
                    react_str += react_res
            else: 
                if "ListOfBonds" in side_list:
                    self.bonds.set_xml(side_list["ListOfBonds"]['Bond'])
                react_res, mol_obj = self.mol_to_str(side_list['ListOfMolecules']['Molecule'])
                mol_obj.set_outer_compt(outer_comp)
                sl.append(mol_obj)
                react_str += react_res
            return react_str, sl
        elif "ProductPattern" in side_xml:
            side_list = side_xml['ProductPattern']
            sl = []
            if '@compartment' in side_list:
                prod_str = "@{}:".format(side_list['@compartment'])
                outer_comp = side_list['@compartment']
            else:
                prod_str = ""
                outer_comp = None
            if isinstance(side_list, list):
                # this is a list of reactants
                for iprod, prod in enumerate(side_list):
                    if "ListOfBonds" in prod:
                        self.bonds.set_xml(prod["ListOfBonds"]['Bond'])
                    if iprod > 0:
                        prod_str += " + "
                    prod_res, mol_obj = self.mol_to_str(prod['ListOfMolecules']['Molecule'])
                    mol_obj.set_outer_compt(outer_comp)
                    sl.append(mol_obj)
                    prod_str += prod_res
            else: 
                sl = []
                if "ListOfBonds" in side_list:
                    self.bonds.set_xml(side_list["ListOfBonds"]['Bond'])
                prod_res, mol_obj = self.mol_to_str(side_list['ListOfMolecules']['Molecule'])
                mol_obj.set_outer_compt(outer_comp)
                sl.append(mol_obj)
                prod_str += prod_res
            return prod_str, sl
        else: 
            print("Can't parse rule XML {}".format(side_xml))

class FuncXML(XMLObj):
    def __init__(self, pattern_xml):
        super().__init__(pattern_xml)

    def resolve_xml(self, func_xml):
        fname = func_xml['@id']
        expression = func_xml['Expression']
        args = []
        if 'ListOfArguments' in func_xml:
            args = self.get_arguments(func_xml['ListOfArguments']['Argument'])
        expr = func_xml['Expression']
        func_str = fname + "("
        if len(args) > 0:
            for iarg, arg in enumerate(args):
                if iarg > 0:
                    func_str += ","
                func_str += arg
        func_str += ")"
        self.item_tuple = (func_str, expression)
        full_str = "{} = {}".format(func_str, expression)
        return full_str 

    def get_arguments(self, arg_xml):
        args = []
        if isinstance(arg_xml, list):
            for arg in arg_xml:
                args.append(arg['@id'])
            return args
        else:
            return [arg_xml['@id']]
###### PATTERNS ###### 
