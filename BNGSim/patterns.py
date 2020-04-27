from BNGSim.bonds import Bonds

###### PATTERNS ###### 
class Pattern:
    def __init__(self, pattern_xml):
        self.molecules = {}
        self.pattern_xml = pattern_xml
        self.bonds = Bonds()
        self.string = self.resolve_xml(self.pattern_xml)

    def __repr__(self):
        return self.string

    def __str__(self):
        return self.string

    def mol_to_str(self, mol_xml):
        # TODO: Document, especially the bond mechanism
        if not isinstance(mol_xml, list):
            # we are going to store molecules, components
            # and compartments in a separate dictionary 
            # for use later
            name = mol_xml['@name'] 
            # we have a single molecule, count 1, order 0
            self.molecules[name] = {'count': 1, 'order': [0]}
            self.molecules[name]['components'] = []
            self.molecules[name]['compartment'] = []
            # start building the string
            mol_str = mol_xml["@name"] + "("
            if "ListOfComponents" in mol_xml:
                # Single molecule can't have bonds
                comp_str = self.comp_to_str(mol_xml["ListOfComponents"]["Component"], mname=mol_xml['@name'])
                mol_str += comp_str
            mol_str += ")"
            if '@compartment' in mol_xml:
                self.molecules[mol_xml['@name']]['compartment'].append(mol_xml['@compartment'])
                mol_str += "@{}".format(mol_xml['@compartment'])
        else:
            # this means we have multiple molecules bonding
            mol_str = ""
            for imol, mol in enumerate(mol_xml):
                name = mol['@name']
                # we can have multiple copies with different 
                # component states bound together. we need to
                # somehow account for it
                if name not in self.molecules:
                    self.molecules[name] = {'count': 1, 'order': [imol]}
                    self.molecules[mol['@name']]['components'] = []
                    self.molecules[mol['@name']]['compartment'] = []
                else:
                    self.molecules[name]['count'] += 1
                    self.molecules[name]['order'].append(imol)
                if imol > 0:
                    # complexing
                    mol_str += "."
                mol_str += mol["@name"] + "("
                if "ListOfComponents" in mol:
                    comp_str = self.comp_to_str(mol['ListOfComponents']['Component'], mname=mol['@name'])
                    mol_str += comp_str
                mol_str += ")"
                if '@compartment' in mol:
                    self.molecules[mol['@name']]['compartment'].append(mol['@compartment'])
                    mol_str += "@{}".format(mol['@compartment'])
        return mol_str

    def comp_to_str(self, comp_xml, mname=None):
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
                if mname is not None:
                    self.molecules[mname]['components'].append(comp_dict)
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
            if mname is not None:
                self.molecules[mname]['components'].append(comp_dict)
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
                if '@compartment' in pattern:
                    obs_str += "@{}:".format(pattern['@compartment'])
                mol = pattern['ListOfMolecules']['Molecule']
                obs_res = self.mol_to_str(mol)
                obs_str += obs_res
        else:
            if "ListOfBonds" in patterns:
                self.bonds.set_xml(patterns["ListOfBonds"]["Bond"])
            mol = patterns['ListOfMolecules']["Molecule"]
            obs_str = ""
            if '@compartment' in patterns:
                obs_str += "@{}:".format(patterns['@compartment'])
            obs_str += self.mol_to_str(mol)
        return obs_str

class SpeciesPattern(Pattern):
    def __init__(self, pattern_xml):
        super().__init__(pattern_xml)

    def resolve_xml(self, spec_xml):
        pattern = spec_xml['ListOfMolecules']['Molecule']
        # this shouldn't be a list
        if isinstance(pattern, list):
            print("species pattern shouldn't be a list")
        else:
            spec_str = self.mol_to_str(pattern)
        return spec_str

class MolTypePattern(Pattern):
    def __init__(self, pattern_xml):
        super().__init__(pattern_xml)

    def resolve_xml(self, molt_xml):
        molt_str = molt_xml['@id'] + "("
        if 'ListOfComponentTypes' in molt_xml:
            comp_dict = molt_xml['ListOfComponentTypes']['ComponentType']
            if '@id' in comp_dict:
                molt_str += comp_dict['@id']
                if "ListOfAllowedStates" in comp_dict:
                    # we have states
                    al_states = comp_dict['ListOfAllowedStates']['AllowedState']
                    if isinstance(al_states, list):
                        for istate, state in enumerate(al_states):
                            molt_str += "~{}".format(state['@id'])
                    else:
                        molt_str += "~{}".format(al_states['@id'])
            else:
                # multiple components
                for icomp, comp in enumerate(comp_dict):
                    if icomp > 0:
                        molt_str += ","
                    molt_str += comp['@id']
                    if "ListOfAllowedStates" in comp:
                        # we have states
                        al_states = comp['ListOfAllowedStates']['AllowedState']
                        if isinstance(al_states, list):
                            for istate, state in enumerate(al_states):
                                molt_str += "~{}".format(state['@id'])
                        else:
                            molt_str += "~{}".format(al_states['@id'])
        molt_str += ")"
        return molt_str

class RulePattern(Pattern):
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
        if 'RateLaw' not in pattern_xml:
            print("Rule seems to be missing a rate law, please make sure that XML exporter of BNGL supports whatever you are doing!")
        rate_law = self.resolve_ratelaw(pattern_xml['RateLaw'])
        # We need to set self.item_tuple to 
        # (rule_name, LHS, RHS, rate_law)
        self.item_tuple = [rule_name, lhs, "->", rhs, rate_law]
        return "{}: {} -> {} {}".format(rule_name, lhs, rhs, rate_law)

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
            import IPython
            IPython.embed()
        return rate_cts

    def resolve_rxn_side(self, side_xml):
        # this is either reactant or product
        if side_xml is None:
            return "0"
        elif 'ReactantPattern' in side_xml:
            # this is a lhs/reactant side
            side_list = side_xml['ReactantPattern']
            if isinstance(side_list, list):
                # this is a list of reactants
                if '@compartment' in side_list:
                    react_str = "@{}:".format(side_list['@compartment'])
                else:
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
                if '@compartment' in side_list:
                    react_str = "@{}:".format(side_list['@compartment'])
                else:
                    react_str = ""
                if "ListOfBonds" in side_list:
                    self.bonds.set_xml(side_list["ListOfBonds"]['Bond'])
                react_res = self.mol_to_str(side_list['ListOfMolecules']['Molecule'])
                react_str += react_res
            return react_str
        elif "ProductPattern" in side_xml:
            side_list = side_xml['ProductPattern']
            if isinstance(side_list, list):
                # this is a list of reactants
                if '@compartment' in side_list:
                    prod_str = "@{}:".format(side_list['@compartment'])
                else:
                    prod_str = ""
                for iprod, prod in enumerate(side_list):
                    if "ListOfBonds" in prod:
                        self.bonds.set_xml(prod["ListOfBonds"]['Bond'])
                    if iprod > 0:
                        prod_str += " + "
                    prod_res = self.mol_to_str(prod['ListOfMolecules']['Molecule'])
                    prod_str += prod_res
            else: 
                if '@compartment' in side_list:
                    prod_str = "@{}:".format(side_list['@compartment'])
                else:
                    prod_str = ""
                if "ListOfBonds" in side_list:
                    self.bonds.set_xml(side_list["ListOfBonds"]['Bond'])
                prod_res = self.mol_to_str(side_list['ListOfMolecules']['Molecule'])
                prod_str += prod_res
            return prod_str
        else: 
            print("Can't parse rule XML {}".format(side_xml))

class FuncPattern(Pattern):
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
