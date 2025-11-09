# spiral chain
from .quadrupole_utils import *


class ChiralAxialType2(ChiralBase):
    def __init__(self, mol):
        super().__init__(mol)

    # whether a bond and a ring share an atom
    def bond_ring(self, bond_, ring_):
        if bond_.GetBeginAtomIdx() in ring_:
            return True
        elif bond_.GetEndAtomIdx() in ring_:
            return True
        else:
            return False

    # find cumulated ene
    def cum_ene(self, double_bonds_):
        # find all bonds that share an atom
        bonds_unit = []
        for idx_1 in range(len(double_bonds_)):
            idx_1_nei = [double_bonds_[idx_1]]
            for idx_2 in range(idx_1 + 1, len(double_bonds_)):
                if self.pub_atom(double_bonds_[idx_1], double_bonds_[idx_2]):
                    idx_1_nei.append(double_bonds_[idx_2])
            bonds_unit.append(idx_1_nei)

        # merge
        bonds_chain_1 = [tuple(chain) for chain in bonds_unit]
        flag = True
        whiles_cnt = [0, 0]
        while flag:
            whiles_cnt[0] += 1
            if whiles_cnt[0] > 200:
                warnings.warn("v2 1st more than 200 times")
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

            # needed, rdkit object to idx
            if self.check_chain_equal(bonds_chain_2) == self.check_chain_equal(bonds_chain_1):
                flag = False
            else:
                # deepcopy fail for bonds chain in rdkit, copy is enough
                bonds_chain_1 = bonds_chain_2.copy()

        flag = True
        bonds_chain_1 = bonds_chain_1[::-1]
        while flag:
            whiles_cnt[1] += 1
            if whiles_cnt[1] > 200:
                warnings.warn("v2 2nd more than 200 times")
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

    # find all spiral chain
    def get_spi_chain(self, spi_):
        # find all spiral atoms that share a ring
        chain_unit = []
        for i in range(len(spi_)):
            i_nei = [spi_[i]]
            for j in range(i + 1, len(spi_)):
                if self.pub_ring(spi_[i], spi_[j]):
                    i_nei.append(spi_[j])
            chain_unit.append(i_nei)

        # merge
        chain_1 = [tuple(i) for i in chain_unit]
        flag = True
        whiles_cnt = [0, 0]
        while flag:
            whiles_cnt[0] += 1
            if whiles_cnt[0] > 200:
                warnings.warn("v2 3rd more than 200 times")
                return []
            chain_2 = []
            for i in range(len(chain_1)):
                k = 0
                for j in range(i + 1, len(chain_1)):  # whether in the same ring
                    l_1 = len(set(chain_1[i]))
                    l_2 = len(set(chain_1[j]))
                    l_3 = len(set(chain_1[i] + chain_1[j]))
                    if l_3 < (l_1 + l_2):
                        delta = set(chain_1[i] + chain_1[j])
                        chain_2.append(tuple(delta))
                    else:
                        k = k + 1
                if k == len(chain_1) - i - 1:
                    chain_2.append(chain_1[i])

            chain_2 = list(set(chain_2))

            if chain_2 == chain_1:
                flag = False
            else:
                chain_1 = chain_2.copy()

        flag = True
        chain_1 = chain_1[::-1]
        while flag:
            whiles_cnt[1] += 1
            if whiles_cnt[1] > 200:
                warnings.warn("v2 4th more than 200 times")
                return []
            chain_2 = []
            for i in range(len(chain_1)):
                k = 0
                for j in range(i + 1, len(chain_1)):
                    l_1 = len(set(chain_1[i]))
                    l_2 = len(set(chain_1[j]))
                    l_3 = len(set(chain_1[i] + chain_1[j]))
                    if l_3 < (l_1 + l_2):
                        delta = set(chain_1[i] + chain_1[j])
                        chain_2.append(tuple(delta))
                    else:
                        k = k + 1
                if k == len(chain_1) - i - 1:
                    chain_2.append(chain_1[i])

            chain_2 = list(set(chain_2))

            if chain_2 == chain_1:
                flag = False
            else:
                chain_1 = chain_2.copy()

        chain_2 = [sorted(i) for i in chain_2]
        return chain_2

    # find end atoms
    def find_end_atom(self, spi_chain_, bonds_chain_2_):
        ssr = Chem.GetSymmSSSR(self.mol)
        end_atoms = []

        for chain in spi_chain_:  # for each chain
            for n in range(1, len(chain)+1):  # for the len_
                for i in range(len(chain) - n + 1):
                    chain_temt = chain[i:i + n + 1]  # get the defined len_
                    chain_atom = []

                    for j in chain_temt:
                        for k in ssr:
                            if j in list(k):
                                chain_atom.extend(list(k))

                    chain_atom = list(set(chain_atom))  # get all atoms about the spiral chain

                    for bonds_chain in bonds_chain_2_:
                        # print(bonds_chain)
                        chain_1 = [bond.GetBeginAtomIdx() for bond in bonds_chain]
                        chain_2 = [bond.GetEndAtomIdx() for bond in bonds_chain]
                        bonds_chain_atom = list(set(chain_1 + chain_2))
                        for atom in bonds_chain_atom:
                            if atom in chain_atom:
                                chain_atom = chain_atom + bonds_chain_atom

                    # print(chain_atom)
                    chain_boundary = []

                    for idx in chain_atom:  # all atoms in a spiral chain
                        atom = self.atoms[idx]
                        atom_nei = atom.GetNeighbors()
                        atom_out_nei = []
                        for nei in atom_nei:
                            if nei.GetIdx() not in chain_atom:
                                atom_out_nei.append(nei.GetIdx())  # outside neighbors of the end atoms
                        # the same or not
                        if len(atom_out_nei) == 2:
                            if self.CIP_list[atom_out_nei[0]] > self.CIP_list[atom_out_nei[1]]:
                                end_atom = [idx, atom_out_nei]  # one end atom and its neighbors
                                chain_boundary.append(end_atom)
                            elif self.CIP_list[atom_out_nei[0]] < self.CIP_list[atom_out_nei[1]]:
                                atom_out_nei = atom_out_nei[::-1]
                                end_atom = [idx, atom_out_nei]  # another end atom and its neighbors
                                chain_boundary.append(end_atom)
                    if len(chain_boundary) >= 2:
                        end_atoms.append([chain_temt, chain_boundary])
        '''
        end_atoms= [
                    [[spiral atom1],[end1,[higher rank neighbor,lower xxx]],[end2,[higher rank neighbor,lower xxx]]...]
                    [[spiral atom2],[end1,[higher rank neighbor,lower xxx]],[end2,[higher rank neighbor,lower xxx]]...]
                    ...
                   ]
        '''
        return end_atoms

    # get the chirality matrix
    def get_chi_mat(self):
        spi_ = self.find_spiral_atoms()
        chain_ = self.get_spi_chain(spi_)
        double_bonds_ = self.get_double_bond()
        bonds_chain_2_ = self.cum_ene(double_bonds_)
        end_atoms_ = self.find_end_atom(chain_, bonds_chain_2_)
        chi_spi_, ends_ = [], []
        mats, dets, norm_cp, signs = [], [], [], []  # for each conf
        neigh_ids = []

        for spi_chain in end_atoms_:
            spi_atoms = spi_chain[0]
            ends = spi_chain[1]
            for i in range(len(ends)):
                boundary = int(ends[i][0])
                nei_1 = ends[i][1]
                for j in range(i + 1, len(ends)):
                    chi_spi_.append(spi_atoms)
                    ends_.append((ends[i], ends[j]))
                    mat_confs = []
                    det_confs = []
                    norm_det_confs = []
                    sign_confs = []
                    for conf_ in self.coordinates:
                        nei_2 = ends[j][1]
                        cri_atoms = [boundary, nei_1[0], nei_1[1], nei_2[0], nei_2[1]]
                        cri_atoms_cor = [conf_[k] for k in cri_atoms]
                        a = cri_atoms_cor[1] - cri_atoms_cor[0]
                        b = cri_atoms_cor[2] - cri_atoms_cor[0]
                        c = cri_atoms_cor[3] - cri_atoms_cor[4]
                        cp_max = np.linalg.norm(np.cross(a, b)) * np.linalg.norm(c)
                        mat = np.array([a, b, c])
                        mat_confs.append(mat)
                        det_, sign_ = self.criterion(mat)
                        det_confs.append(det_)
                        norm_det_confs.append(det_/cp_max)
                        sign_confs.append(sign_)
                    mats.append(mat_confs)
                    dets.append(det_confs)
                    norm_cp.append(norm_det_confs)
                    signs.append(sign_confs)
                    neigh_ids.append(cri_atoms[1:])

        axes_label = []
        for i in range(len(chi_spi_)):
            axes_label.extend(chi_spi_[i])
            axes_label.append(ends_[i][0][0])
            axes_label.append(ends_[i][1][0])
        """axes atoms, more than other results!"""
        return {"spiral id": chi_spi_, "ends": ends_, "chiral axes": [(one,) for one in set(axes_label)], 
                "quadrupole matrix": mats, "determinant": dets, "norm CP": norm_cp, "sign": signs, "neighbor ids": neigh_ids}
