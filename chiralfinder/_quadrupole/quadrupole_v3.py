# Allene-like
from .quadrupole_utils import *


class ChiralAxialType3(ChiralBase):
    def __init__(self, mol):
        super().__init__(mol)

    # find cumulated ene
    def cum_ene(self, double_bonds_):
        # find all bonds that share an atom
        bonds_unit = []
        for idx_1 in range(len(double_bonds_)):
            idx_1_nei = set()
            idx_1_nei.add(double_bonds_[idx_1])
            for idx_2 in range(idx_1 + 1, len(double_bonds_)):
                if self.pub_atom(double_bonds_[idx_1], double_bonds_[idx_2]):
                    idx_1_nei.add(double_bonds_[idx_2])
            bonds_unit.append(list(idx_1_nei))

        # bonds_chain_1 = bonds_unit
        bonds_chain_1 = [tuple(chain) for chain in bonds_unit]
        bonds_chain_1 = list(set(bonds_chain_1))

        # more than one double bond
        bonds_chain_3 = []
        for bonds_chain in bonds_chain_1:
            if len(bonds_chain) != 1:
                bonds_chain_3.append(bonds_chain)
        bonds_chain_1 = bonds_chain_3.copy()

        # merge
        flag = True
        whiles_cnt = [0, 0]
        while flag:
            whiles_cnt[0] += 1
            if whiles_cnt[0] > 200:
                warnings.warn("v3 1st more than 200 times")
                return []
            bonds_chain_2 = []
            for i in range(len(bonds_chain_1)):
                k = 0
                for j in range(i + 1, len(bonds_chain_1)):
                    l_1 = len(set(bonds_chain_1[i]))
                    l_2 = len(set(bonds_chain_1[j]))
                    l_3 = len(set(bonds_chain_1[i] + bonds_chain_1[j]))
                    if l_3 < (l_1 + l_2):
                        delta = set(bonds_chain_1[i] + bonds_chain_1[j])
                        bonds_chain_2.append(tuple(delta))
                    else:
                        k = k + 1
                if k == len(bonds_chain_1) - i - 1:
                    bonds_chain_2.append(bonds_chain_1[i])

            bonds_chain_2 = list(set(bonds_chain_2))
            # print(bonds_chain_2)
            if self.check_chain_equal(bonds_chain_2) == self.check_chain_equal(bonds_chain_1):
                flag = False
            else:
                bonds_chain_1 = bonds_chain_2.copy()

        flag = True
        bonds_chain_1 = bonds_chain_1[::-1]
        while flag:
            whiles_cnt[1] += 1
            if whiles_cnt[1] > 200:
                warnings.warn("v3 2nd more than 200 times")
                return []
            bonds_chain_2 = []
            for i in range(len(bonds_chain_1)):
                k = 0
                for j in range(i + 1, len(bonds_chain_1)):
                    l_1 = len(set(bonds_chain_1[i]))
                    l_2 = len(set(bonds_chain_1[j]))
                    l_3 = len(set(bonds_chain_1[i] + bonds_chain_1[j]))
                    if l_3 < (l_1 + l_2):
                        delta = set(bonds_chain_1[i] + bonds_chain_1[j])
                        bonds_chain_2.append(tuple(delta))
                    else:
                        k = k + 1
                if k == len(bonds_chain_1) - i - 1:
                    bonds_chain_2.append(bonds_chain_1[i])

            bonds_chain_2 = list(set(bonds_chain_2))

            if self.check_chain_equal(bonds_chain_2) == self.check_chain_equal(bonds_chain_1):
                flag = False
            else:
                bonds_chain_1 = bonds_chain_2.copy()

        return bonds_chain_2

    # find the end atoms
    def find_ene_end_atoms(self, cum_ene_):
        ssr = Chem.GetSymmSSSR(self.mol)

        cum_ene_end = []
        for ene in cum_ene_:
            ene_atom_idx_1 = [bond.GetBeginAtomIdx() for bond in ene]
            ene_atom_idx_2 = [bond.GetEndAtomIdx() for bond in ene]
            ene_atom_idx = list(set(ene_atom_idx_1 + ene_atom_idx_2))  # get all atoms
            ene_atom_idx = sorted(ene_atom_idx)
            ends_idx = [ene_atom_idx[0], ene_atom_idx[-1]]  # get the end atoms

            # whether with another ring
            for idx in ends_idx:
                for ring in ssr:
                    if idx in list(ring):
                        ene_atom_idx = list(set(ene_atom_idx + list(ring)))
                        ends_idx = list(set(ends_idx + list(ring)))

            # outside neighbors
            ene_end = []
            for end_atom_idx in ends_idx:
                end_atom = self.atoms[end_atom_idx]
                end_atom_nei = end_atom.GetNeighbors()
                end_out_nei = []
                for nei in end_atom_nei:
                    if nei.GetIdx() not in ene_atom_idx:
                        end_out_nei.append(nei.GetIdx())

                # whether two different atoms
                if len(end_out_nei) == 2:
                    if self.CIP_list[end_out_nei[0]] > self.CIP_list[end_out_nei[1]]:
                        ene_end.append([end_atom_idx, end_out_nei])
                    elif self.CIP_list[end_out_nei[0]] < self.CIP_list[end_out_nei[1]]:
                        ene_end.append([end_atom_idx, end_out_nei[::-1]])

                if len(ene_end) == 2:
                    if ene_end not in cum_ene_end:
                        cum_ene_end.append(ene_end)

        '''
        cum_ene_end=[
                    [end1,[higher rank neighbor,lower xxx],end2,[higher rank neighbor,lower xxx]]
                    [end1,[higher rank neighbor,lower xxx],end2,[higher rank neighbor,lower xxx]]
                    ...
                    ]
        '''
        return cum_ene_end

    # get
    def get_chi_mat(self):
        double_bonds_ = self.get_double_bond()
        cum_ene_ = self.cum_ene(double_bonds_)
        cum_ene_end_ = self.find_ene_end_atoms(cum_ene_)

        chi_axial_ = []
        mats, dets, norm_cp, signs = [], [], [], []  # for each conf
        neigh_ids = []
        for ene in cum_ene_end_:
            end_1 = ene[0]
            end_2 = ene[1]
            chi_axial_.append((end_1[0], end_2[0]))
            mat_confs = []
            det_confs = []
            norm_det_confs = []
            sign_confs = []
            for conf_ in self.coordinates:
                end_cor_1 = [conf_[end_1[0]]]
                end_cor_2 = [conf_[i] for i in end_1[1]]
                end_cor_3 = [conf_[i] for i in end_2[1]]
                end_cor = end_cor_1 + end_cor_2 + end_cor_3

                a = end_cor[1] - end_cor[0]
                b = end_cor[2] - end_cor[0]
                c = end_cor[4] - end_cor[3]
                cp_max = np.linalg.norm(np.cross(a, b)) * np.linalg.norm(c)
                mat = np.array([a, b, c])
                mat_confs.append(mat)
                det_, sign_ = self.criterion(mat)
                det_confs.append(det_)
                norm_det_confs.append(det_/cp_max)
                sign_confs.append(sign_)
                neigh_ids.append(end_1[1] + end_2[1])
            mats.append(mat_confs)
            dets.append(det_confs)
            norm_cp.append(norm_det_confs)
            signs.append(sign_confs)
        return {"axial id": chi_axial_, "chiral axes": chi_axial_, "quadrupole matrix": mats, 
                "determinant": dets, "norm CP": norm_cp, "sign": signs, "neighbor ids": neigh_ids}
