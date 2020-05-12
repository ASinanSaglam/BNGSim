from BNGSim.patterns import ObsPattern, MolTypePattern, RulePattern, FuncPattern, SpeciesPattern

###### MODEL STRUCTURES ###### 
# Objects in the model
class ModelBlock:
    def __init__(self):
        self._item_dict = {}

    def __len__(self):
        return len(self._item_dict)

    def __repr__(self):
        # overwrites what the class representation
        # shows the items in the model block in 
        # say ipython
        return str(self._item_dict)

    def __getitem__(self, key):
        if isinstance(key, int):
            # get the item in order
            return list(self._item_dict.keys())[key]
        return self._item_dict[key]

    def __setitem__(self, key, value):
        self._item_dict[key] = value

    def __delitem__(self, key):
        if key in self._item_dict:
            self._item_dict.pop(key)
        else: 
            print("Item {} not found".format(key))

    def __iter__(self):
        return self._item_dict.keys().__iter__()

    def __contains__(self, key):
        return key in self._item_dict

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


# TODO: Add a LOT of error handling
class Parameters(ModelBlock):
    '''
    Class containing parameters
    '''
    def __init__(self):
        self.expressions = {}
        self.values = {}
        super().__init__()
        self.name = "parameters"

    def __setattr__(self, name, value):
        changed = False
        if hasattr(self, "_item_dict"):
            if name in self._item_dict.keys():
                try: 
                    new_value = float(value)
                    changed = True
                    self._item_dict[name] = new_value
                except:
                    self._item_dict[name] = value
        if changed:
            self.__dict__[name] = new_value
        else:
            self.__dict__[name] = value

    def __str__(self):
        # overwrites what the method returns when 
        # it's converted to string
        block_lines = ["\nbegin {}".format(self.name)]
        for item in self._item_dict.keys():
            block_lines.append("  " + "{} {}".format(item, self._item_dict[item]))
        block_lines.append("end {}\n".format(self.name))
        return "\n".join(block_lines)

    def add_item(self, item_tpl):
        name, value = item_tpl
        self._item_dict[name] = value
        try:
            setattr(self, name, value)
        except:
            print("can't set {} to {}".format(name, value))
            pass

    def parse_xml_block(self, block_xml):
        # 
        if isinstance(block_xml, list):
            for b in block_xml:
                self.values[b['@id']] = b['@value']
                if '@expr' in b:
                    self.expressions[b['@id']] = b['@expr']
                    self.add_item((b['@id'],b['@expr']))
                else:
                    self.add_item((b['@id'],b['@value']))
        else:
            self.values[block_xml['@id']] = block_xml['@value']
            if '@expr' in block_xml:
                self.expressions[block_xml['@id']] = block_xml['@expr']
                self.add_item((block_xml['@id'], block_xml['@expr']))
            else:
                self.add_item((block_xml['@id'], block_xml['@value']))
        # 

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

    def __getitem__(self, key):
        if isinstance(key, str):
            # our keys are pattern objects
            for ikey in self._item_dict:
                if key == ikey.string:
                    return self._item_dict[ikey]
        if isinstance(key, int):
            # get the item in order
            return list(self._item_dict.keys())[key]
        return self._item_dict[key]

    def __setitem__(self, key, value):
        for ikey in self._item_dict:
            if key == ikey.string:
                self._item_dict[ikey] = value

    def __contains__(self, key):
        for ikey in self._item_dict:
            if key == ikey.string:
                return True
        return False

    def parse_xml_block(self, block_xml):
        #TODO: Eventually regenerate patterns 
        # with bond handling from XML instead of 
        # reading the name directly to stay consistent
        #TODO This is especially important since XML dumps 
        # compartments with format @X::Species and that 
        # seems uncommon currently
        if isinstance(block_xml, list):
            for sd in block_xml:
                pattern = SpeciesPattern(sd)
                self.add_item((pattern,sd['@concentration']))
        else:
            pattern = SpeciesPattern(block_xml)
            self.add_item((pattern,block_xml['@concentration']))

    def add_item(self, item_tpl):
        name, val = item_tpl
        self._item_dict[name] = val

class MoleculeTypes(ModelBlock):
    '''
    Class containing molecule types 
    '''
    def __init__(self):
        super().__init__()
        self.name = "molecule types"

    def __repr__(self):
        return str(list(self._item_dict.keys()))

    def add_item(self, item_tpl):
        name, = item_tpl
        self._item_dict[name] = ""

    def __getitem__(self, key):
        if isinstance(key, str):
            # our keys are pattern objects
            for ikey in self._item_dict:
                if key == ikey.string:
                    return self._item_dict[ikey]
        if isinstance(key, int):
            # get the item in order
            return list(self._item_dict.keys())[key]
        return self._item_dict[key]

    def __setitem__(self, key, value):
        for ikey in self._item_dict:
            if key == ikey.string:
                self._item_dict[ikey] = value

    def __contains__(self, key):
        for ikey in self._item_dict:
            if key == ikey.string:
                return True
        return False

    def __str__(self):
        # overwrites what the method returns when 
        # it's converted to string
        block_lines = ["\nbegin {}".format(self.name)]
        for item in self._item_dict.keys():
            block_lines.append("  " + "{}".format(item))
        block_lines.append("end {}\n".format(self.name))
        return "\n".join(block_lines)

    def parse_xml_block(self, block_xml):
        if isinstance(block_xml, list):
            for md in block_xml:
                pattern = MolTypePattern(md)
                self.add_item((pattern,))
        else:
            pattern = MolTypePattern(block_xml)
            self.add_item((pattern,))


class Observables(ModelBlock):
    # TODO: Compartments
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

    def __getitem__(self, key):
        if isinstance(key, int):
            # get the item in order
            return self._item_dict[list(self._item_dict.keys())[key]][1]
        return self._item_dict[key]

    def parse_xml_block(self, block_xml):
        #
        if isinstance(block_xml, list):
            for b in block_xml:
                pattern = ObsPattern(b['ListOfPatterns'])
                self.add_item((b['@type'], b['@name'], pattern))
        else: 
            pattern = ObsPattern(block_xml['ListOfPatterns'])
            self.add_item((block_xml['@type'], block_xml['@name'], pattern))
        # 


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
        for item in self._item_dict.keys():
            block_lines.append("  " + 
                    "{} = {}".format(item, self._item_dict[item]))
        block_lines.append("end {}\n".format(self.name))
        return "\n".join(block_lines)

    def parse_xml_block(self, block_xml):
        if isinstance(block_xml, list):
             for func in block_xml:
                 fpatt = FuncPattern(func)
                 self.add_item(fpatt.item_tuple)
        else:
             fpatt = FuncPattern(block_xml)
             self.add_item(fpatt.item_tuple)

class Compartments(ModelBlock):
    '''
    Class for compartments
    '''
    def __init__(self):
        super().__init__()
        self.name = "compartments"

    def __str__(self):
        # overwrites what the method returns when 
        # it's converted to string
        block_lines = ["\nbegin {}".format(self.name)]
        for item in self._item_dict.keys():
            comp_line = "  {} {} {}".format(item, 
                            self._item_dict[item][0],
                            self._item_dict[item][1])
            if self._item_dict[item][2] is not None:
                comp_line += " {}".format(self._item_dict[item][2])
            block_lines.append(comp_line)
        block_lines.append("end {}\n".format(self.name))
        return "\n".join(block_lines)

    def add_item(self, item_tpl):
        name, dim, size, outside = item_tpl
        self._item_dict[name] = [dim, size, outside]

    def parse_xml_block(self, block_xml):
        # 
        if isinstance(block_xml, list):
            for comp in block_xml:
                cname = comp['@id']
                dim = comp['@spatialDimensions']
                size = comp['@size']
                if '@outside' in comp:
                    outside = comp['@outside']
                else:
                    outside = None
                self.add_item( (cname, dim, size, outside) )
        else:
            cname = block_xml['@id']
            dim = block_xml['@spatialDimensions']
            size = block_xml['@size']
            if '@outside' in block_xml:
                outside = block_xml['@outside']
            else:
                outside = None
            self.add_item( (cname, dim, size, outside) )
        #

class Rules(ModelBlock):
    def __init__(self):
        super().__init__()
        self.name = "reaction rules"

    def __str__(self):
        # TODO: printing also needs a lot of adjusting
        block_lines = ["\nbegin {}".format(self.name)]
        for item in self._item_dict.keys():
            block_lines.append(str(self._item_dict[item]))
        block_lines.append("end {}\n".format(self.name))
        return "\n".join(block_lines)

    def parse_xml_block(self, block_xml):
        if isinstance(block_xml, list):
            for rd in block_xml:
                rpattern = RulePattern(rd)
                self.add_item((rpattern.name, rpattern))
        else:
            rpattern = RulePattern(block_xml)
            self.add_item((rpattern.name, rpattern))
        self.consolidate_rules()

    def consolidate_rules(self):
        '''
        Generated XML only has unidirectional rules
        and uses "_reverse_" tag to make bidirectional 
        rules for NFSim. Take all the "_reverse_" tagged
        rules and convert them to bidirectional rules
        '''
        delete_list = []
        for item_key in self._item_dict:
            rxn_pat = self._item_dict[item_key]
            if item_key.startswith("_reverse_"):
                # this is the reverse of another reaction
                reverse_of = item_key.replace("_reverse_", "")
                # ensure we have the original
                if reverse_of in self._item_dict:
                    # make bidirectional and add rate law
                    r1 = self._item_dict[reverse_of].rate_law[0]
                    r2 = rxn_pat.rate_law[0]
                    self._item_dict[reverse_of].set_rate_law((r1,r2))
                    # mark reverse for deletion
                    delete_list.append(item_key)
        # delete items marked for deletion
        for del_item in delete_list:
            self._item_dict.pop(del_item)

class Actions(ModelBlock):
    def __init__(self):
        super().__init__()
        self.name = "actions"
        self._action_list = ["generate_network", "generate_hybrid_model","simulate", "simulate_ode", "simulate_ssa", "simulate_pla", "simulate_nf", "parameter_scan", "bifurcate", "readFile", "writeFile", "writeModel", "writeNetwork", "writeXML", "writeSBML", "writeMfile", "writeMexfile", "writeMDL", "visualize", "setConcentration", "addConcentration", "saveConcentration", "resetConcentrations", "setParameter", "saveParameters", "resetParameters", "quit", "setModelName", "substanceUnits", "version", "setOption"]

    def __str__(self):
        # TODO: figure out every argument that has special 
        # requirements, e.g. method requires the value to 
        # be a string
        block_lines = []
        for item in self._item_dict.keys():
            action_str = "{}(".format(item) + "{"
            for iarg,arg in enumerate(self._item_dict[item]):
                val = arg[1]
                arg = arg[0]
                if iarg > 0:
                    action_str += ","
                if arg == "method":
                    action_str += '{}=>"{}"'.format(arg, val)
                else:
                    action_str += '{}=>{}'.format(arg, val)
            action_str += "})"
            block_lines.append(action_str)
        return "\n".join(block_lines)

    def add_action(self, action_type, action_args):
        '''
        adds action, needs type as string and args as list of tuples
        (which preserve order) of (argument, value) pairs
        '''
        if action_type in self._action_list:
            self._item_dict[action_type] = action_args
        else:
            print("Action type {} not valid".format(action_type))
###### MODEL STRUCTURES ###### 
