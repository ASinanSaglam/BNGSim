class Molecule:
    def __init__(self,  mol_list):
        self.mol_list = mol_list

    def __str__(self):
        pat_str = ""
        for imol, mol in enumerate(self.mol_list):
            if imol > 0:
                pat_str += "."
            mol_str = mol["name"]
            if imol == 0:
                if mol['outer_comp'] is not None:
                    mol_str = "@{}:{}".format(mol['outer_comp'],mol_str)
            if len(mol["components"]) > 0:
                mol_str += "("
                for icomp, comp in enumerate(mol["components"]):
                    if icomp > 0:
                        mol_str += ","
                    comp_str = comp["name"]
                    # only for moltypes
                    if "states" in comp:
                        for istate, state in enumerate(comp["states"]):
                            comp_str += "~{}".format(state)
                    # for any other pattern
                    if "state" in comp:
                        comp_str += "~{}".format(comp["state"])
                    if "bonds" in comp:
                        for bond in comp["bonds"]:
                            comp_str += "!{}".format(bond)
                    mol_str += comp_str 
                mol_str += ")"
            pat_str += mol_str
            if mol["compartment"] is not None:
                pat_str += "@{}".format(mol["compartment"])
        return pat_str

    def set_outer_compt(self, out_comp):
        self.mol_list[0]['outer_comp'] = out_comp

    def add_component(self, name, state=None, states=None, num=None):
        for imol, mol in enumerate(self.mol_list):
            if num is not None:
                if imol == num:
                    comp_dict = {"name": name}
                    if state is not None:
                        comp_dict["state"] = state
                    if states is not None:
                        comp_dict["states"] = states
                    mol['components'].append(comp_dict)
            else:
                comp_dict = {"name": name}
                if state is not None:
                    comp_dict["state"] = state
                if states is not None:
                    comp_dict["states"] = states
                mol['components'].append(comp_dict)

    def set_compartment(self, compt, num=None):
        for imol, mol in enumerate(self.mol_list):
            if num is not None:
                if imol == num:
                    mol['compartment'] = compt
            else:
                mol['compartment'] = compt
