
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
        # self.bonds_dict is a dictionary you can key
        # with the tuple taken from the ID and then 
        # get a bond ID cleanly
        if isinstance(bonds_xml, list):
            for ibond, bond in enumerate(bonds_xml): 
                bond_partner_1, bond_partner_2 = self.tpls_from_bond(bond)
                if bond_partner_1 not in self.bonds_dict:
                    self.bonds_dict[bond_partner_1] = [ibond+1]
                else:
                    self.bonds_dict[bond_partner_1].append([ibond+1])
                if bond_partner_2 not in self.bonds_dict:
                    self.bonds_dict[bond_partner_2] = [ibond+1]
                else:
                    self.bonds_dict[bond_partner_2].append(ibond+1)
        else:
            bond_partner_1, bond_partner_2 = self.tpls_from_bond(bonds_xml)
            self.bonds_dict[bond_partner_1] = [1]
            self.bonds_dict[bond_partner_2] = [1]
###### BONDS #####

