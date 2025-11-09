# Quadrupole Matrix, spiral atom
from .quadrupole_utils import *


class ChiralAxialType1(ChiralBase):
    """axial spiral (1 atom)"""
    def __init__(self, mol):
        super().__init__(mol)

    # find spiral atoms with "axial chirality"
    def find_chi_spi(self, spi_):
        chi_spi = []
        ssr = Chem.GetSymmSSSR(self.mol)
        for i in spi_:
            # get the rings where spiral atoms are in
            i_ring = []
            for ring in ssr:
                if i in ring:
                    i_ring.append(list(ring))
            # classify the spiral atom's neighbors with rings
            nei_1 = self.atoms[i].GetNeighbors()
            nei_1 = [x.GetIdx() for x in nei_1]
            raw_nei_1 = deepcopy(nei_1)
            nei_2 = []
            for j in range(len(raw_nei_1)):
                if raw_nei_1[j] in i_ring[0]:
                    nei_2.append(raw_nei_1[j])
                    nei_1.remove(raw_nei_1[j])
            # whether neighbors who are in the same ring are equivalent
            if (nei_1[0] - nei_1[1]) * (nei_2[0] - nei_2[1]) != 0:
                chi_spi.append([i, nei_1, nei_2])

        return chi_spi

    # get the chiral matrices
    def get_chi_mat(self):
        chi_spi_ = self.find_chi_spi(self.find_spiral_atoms())

        mats, dets, norm_cp, signs = [], [], [], []  # for each conf
        neigh_ids = []
        for one in chi_spi_:
            # (id_, rank_), sort by rank, increasing
            nei_1_id_rank = [(i, self.CIP_list[self.atoms[i].GetIdx()]) for i in one[1]]
            nei_2_id_rank = [(i, self.CIP_list[self.atoms[i].GetIdx()]) for i in one[2]]
            nei_1_id_rank = sorted(nei_1_id_rank, key=lambda x: x[1])
            nei_2_id_rank = sorted(nei_2_id_rank, key=lambda x: x[1])
            neigh_id_rank = nei_1_id_rank + nei_2_id_rank

            mat_confs = []
            det_confs = []
            norm_det_confs = []
            sign_confs = []
            for conf_ in self.coordinates:
                neigh_cor = [conf_[one[0]]]
                for one_ in neigh_id_rank:
                    neigh_cor.append(conf_[one_[0]])
                # get the matrix
                a = neigh_cor[1] - neigh_cor[0]
                b = neigh_cor[2] - neigh_cor[0]
                c = neigh_cor[3] - neigh_cor[4]
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
            neigh_ids.append([x[0] for x in neigh_id_rank])
        return {"spiral id": chi_spi_, "chiral axes": [(one[0],) for one in chi_spi_], "quadrupole matrix": mats, 
                "determinant": dets, "norm CP": norm_cp, "sign": signs, "neighbor ids": neigh_ids}