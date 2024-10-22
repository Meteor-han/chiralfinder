from .quadrupole_utils import *


class ChiralAxialType4(ChiralBase):
    def __init__(self, mol):
        super().__init__(mol)

    def find_connecting_atoms(self, bond_list, atom, prior_atom=None):
        # find all atoms connect to the input 'atom', remove 'prior_atom' in them
        connecting_atoms = []
        for bond in bond_list:
            if bond[0] == atom and not bond[1] == prior_atom:
                connecting_atoms.append(bond[1])
        return connecting_atoms

    def get_chi_mat(self):
        # all bonds
        bond_list = []
        for bond in self.mol.GetBonds():
            a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            # bi-direction
            bond_list.extend([[a1, a2], [a2, a1]])
        # bond connects ph
        phConnecting = []
        for bond in bond_list:
            if self.atoms[bond[0]].GetIsAromatic() and self.atoms[bond[1]].GetIsAromatic():
                if not self.mol.GetBondBetweenAtoms(bond[0], bond[1]).GetIsAromatic():
                    # just remain one direction
                    if bond[0] < bond[1]:
                        phConnecting.append((bond[0], bond[1]))

        # just an order, not CIP, unable to relate to R/S
        res = list(Chem.CanonicalRankAtoms(self.mol, breakTies=False, includeChirality=True, includeIsotopes=True))
        atoms = self.mol.GetAtoms()
        # get the axis
        chiral_axes = []  # merge all confs
        mats, dets, signs = [], [], []  # for each conf
        for bond in phConnecting:
            # phbonds neighbors
            begin_neighbor = [atom.GetIdx() for atom in atoms[bond[0]].GetNeighbors()]
            end_neighbor = [atom.GetIdx() for atom in atoms[bond[1]].GetNeighbors()]
            # outside neighbors
            if bond[1] in begin_neighbor:
                begin_neighbor.remove(bond[1])
            if bond[0] in end_neighbor:
                end_neighbor.remove(bond[0])

            # only consider this type, two atoms (include H)
            if (len(begin_neighbor) != 2) or (len(end_neighbor) != 2):
                continue

            # different then chiral
            if (res[begin_neighbor[0]] != res[begin_neighbor[1]]) and (res[end_neighbor[0]] != res[end_neighbor[1]]):
                chiral_axes.append(bond)

                mat_confs = []
                det_confs = []
                sign_confs = []
                for conf_ in self.coordinates:
                    begin_cor = conf_[bond[0]]
                    end_cor = conf_[bond[1]]

                    # sort the outside neighbors
                    if res[begin_neighbor[0]] > res[begin_neighbor[1]]:
                        begin_1_cor = conf_[begin_neighbor[0]]
                        begin_2_cor = conf_[begin_neighbor[1]]
                    else:
                        begin_1_cor = conf_[begin_neighbor[1]]
                        begin_2_cor = conf_[begin_neighbor[0]]

                    if res[end_neighbor[0]] > res[end_neighbor[1]]:
                        end_1_cor = conf_[end_neighbor[0]]
                        end_2_cor = conf_[end_neighbor[1]]
                    else:
                        end_1_cor = conf_[end_neighbor[1]]
                        end_2_cor = conf_[end_neighbor[0]]

                    a = begin_1_cor - (begin_cor + end_cor) / 2
                    b = begin_2_cor - (begin_cor + end_cor) / 2
                    c = end_1_cor - end_2_cor
                    mat = np.array([a, b, c])

                    mat_confs.append(mat)
                    det_, sign_ = self.criterion(mat)
                    det_confs.append(det_)
                    sign_confs.append(sign_)
                mats.append(mat_confs)
                dets.append(det_confs)
                signs.append(sign_confs)
        return {"chiral axes": chiral_axes, "quadrupole matrix": mats, "determinant": dets, "sign": signs}
