from BNGSim.pattern import Pattern, Molecule, Bonds

###### XMLObjs ###### 
class XMLObj:
    def __init__(self, xml):
        self.items = []
        self.xml = xml
        self.resolve_xml(self.xml)

    def __repr__(self):
        return self.gen_string()

    def __str__(self):
        return self.gen_string()

    def gen_string(self):
        # TODO: Eventually slap in a mechanism where
        # this only runs when the underlying structure 
        # changes 
        sstr = ""
        for iitem, item in enumerate(self.items):
            if iitem > 0:
                sstr += "."
            sstr += str(item)
        return sstr

class ObsXML(XMLObj):
    '''
    An observable is a list of patterns where a pattern
    is a list of molecules
    '''
    def __init__(self, xml):
        super().__init__(xml)

    def gen_string(self):
        obs_str = ""
        for iobs, obs in enumerate(self.items):
            if iobs > 0:
                obs_str += ","
            obs_str += str(obs)
        return obs_str

    def resolve_xml(self, obs_xml):
        patterns = obs_xml['Pattern']
        if isinstance(patterns, list):
            # we have multiple patterns so this is a list
            for ipattern, pattern in enumerate(patterns): 
                # 
                self.items.append(Pattern(pattern))
        else:
            self.items.append(Pattern(patterns))

class SpeciesXML(XMLObj):
    '''
    A species is a list of molecules
    '''
    def __init__(self, xml):
        super().__init__(xml)

    def __getitem__(self,key):
        return self.pattern[key]

    def gen_string(self):
        return str(self.pattern)

    def resolve_xml(self, spec_xml):
        self.pattern = Pattern(spec_xml)

class MolTypeXML(XMLObj):
    def __init__(self, xml):
        super().__init__(xml)

    def gen_string(self):
        return str(self.molecule)

    def resolve_xml(self, molt_xml):
        molt_name = molt_xml['@id'] 
        mol_dict = {"name": molt_name, "components": [], "compartment": None}
        if 'ListOfComponentTypes' in molt_xml:
            comp_dict = molt_xml['ListOfComponentTypes']['ComponentType']
            if '@id' in comp_dict:
                cd = {"name": comp_dict["@id"]}
                if "ListOfAllowedStates" in comp_dict:
                    # we have states
                    al_states = comp_dict['ListOfAllowedStates']['AllowedState']
                    cd["states"] = []
                    if isinstance(al_states, list):
                        for istate, state in enumerate(al_states):
                            cd["states"].append(state["@id"])
                    else:
                        cd["states"].append(state["@id"])
                mol_dict["components"].append(cd)
            else:
                # multiple components
                for icomp, comp in enumerate(comp_dict):
                    cd = {"name": comp['@id']}
                    if "ListOfAllowedStates" in comp:
                        # we have states
                        al_states = comp['ListOfAllowedStates']['AllowedState']
                        cd['states'] = []
                        if isinstance(al_states, list):
                            for istate, state in enumerate(al_states):
                                cd['states'].append(state['@id'])
                        else:
                            cd['states'].append(state['@id'])
                    mol_dict['components'].append(cd)
        self.molecule = Molecule(mol_dict)

class RuleXML(XMLObj):
    '''
    A rule is a tuple (list of reactant patterns, list of 
    product patterns, list of rate constant functions)
    '''
    def __init__(self, pattern_xml):
        self.bidirectional = False
        super().__init__(pattern_xml)

    def gen_string(self):
        if self.bidirectional:
            return "{}: {} <-> {} {},{}".format(self.name, self.side_string(self.lhs), self.side_string(self.rhs), self.rate_law[0], self.rate_law[1])
        else:
            return "{}: {} -> {} {}".format(self.name, self.side_string(self.lhs), self.side_string(self.rhs), self.rate_law[0])

    def side_string(self, side):
        pat_str = ""
        for ipat, pat in enumerate(side):
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
        self.lhs = self.resolve_rxn_side(pattern_xml['ListOfReactantPatterns'])
        self.rhs = self.resolve_rxn_side(pattern_xml['ListOfProductPatterns'])
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
            return [None]
        elif 'ReactantPattern' in side_xml:
            # this is a lhs/reactant side
            sl = []
            side = side_xml['ReactantPattern']
            # FIXME: how to tackle this?
            if '@compartment' in side:
                outer_comp = side['@compartment']
            else:
                outer_comp = None
            if isinstance(side, list):
                # this is a list of reactants
                for ireact, react in enumerate(side):
                    sl.append(Pattern(react))
            else: 
                sl.append(Pattern(side))
            return sl
        elif "ProductPattern" in side_xml:
            side = side_xml['ProductPattern']
            sl = []
            # FIXME: how to tackle this?
            if '@compartment' in side:
                outer_comp = side['@compartment']
            else:
                outer_comp = None
            if isinstance(side, list):
                # this is a list of reactants
                for iprod, prod in enumerate(side):
                    sl.append(Pattern(prod))
            else: 
                sl.append(Pattern(side))
            return sl
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
