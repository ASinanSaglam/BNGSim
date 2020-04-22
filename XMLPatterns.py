from XMLBonds import Bonds

###### PATTERNS ###### 
class Pattern:
    def __init__(self, pattern_xml):
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

    def resolve_xml(self, molt_xml):
        molt_str = molt_xml['@id'] + "("
        if 'ListOfComponentTypes' in molt_xml:
            comp_dict = molt_xml['ListOfComponentTypes']['ComponentType']
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
        self.item_tuple = [rule_name, lhs, "->", rhs, rate_law]
        return "{}: {} -> {} {}".format(rule_name, lhs, rhs, rate_law)

    def resolve_ratelaw(self, rate_xml):
        rate_type = rate_xml['@type']
        if rate_type == 'Ele':
            rate_cts_xml = rate_xml['ListOfRateConstants']
            rate_cts = rate_cts_xml['RateConstant']['@value']
        elif rate_type == 'Function':
            rate_cts = rate_xml['@name']
        else:
            print("don't recognize rate law type")
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
