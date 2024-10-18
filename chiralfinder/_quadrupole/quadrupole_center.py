from quadrupole_utils import *


class ChiralCenter(ChiralBase):
    def __init__(self, mol, mol_wo_Hs=None, CIP=True):
        super().__init__(mol, mol_wo_Hs, CIP)

    # get the chiral matrices
    def get_chi_mat(self):
        cen_chi = self.find_center_atoms()

        mats, dets, signs = [], [], []  # for each conf
        for i in cen_chi:
            mat_confs = []
            det_confs = []
            sign_confs = []
            for conf_ in self.coordinates:
                atom = self.atoms[i]
                neighbors = atom.GetNeighbors()
                # (id_, rank_), sort by rank, increasing
                neigh_id_rank = [(atom_.GetIdx(), self.CIP_list[atom_.GetIdx()]) for atom_ in neighbors]
                neigh_id_rank = sorted(neigh_id_rank, key=lambda x: x[1])

                # center, atom 1, 2, 3, 4
                neigh_cor = [conf_[i]]
                for one in neigh_id_rank:
                    neigh_cor.append(conf_[one[0]])
                # get the matrix
                if len(neigh_cor) == 4:
                    neigh_cor.insert(1, (neigh_cor[1]+neigh_cor[2]+neigh_cor[3]-neigh_cor[0]*3)/3*-1.0+neigh_cor[0])
                if len(neigh_cor) < 5:
                    continue
                a = neigh_cor[1] - neigh_cor[0]
                b = neigh_cor[4] - neigh_cor[3]
                c = neigh_cor[4] - neigh_cor[2]
                mat = np.array([a, b, c])
                mat_confs.append(mat)
                det_, sign_ = self.criterion(mat)
                det_confs.append(det_)
                sign_confs.append(sign_)
            mats.append(mat_confs)
            dets.append(det_confs)
            signs.append(sign_confs)
        return {"center id": cen_chi, "quadrupole matrix": mats, "determinant": dets, "sign": signs}
