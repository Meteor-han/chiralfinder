from .quadrupole_utils import *
from .quadrupole_center import ChiralCenter


class ChiralAxialType6(ChiralCenter):
    def __init__(self, mol, mol_wo_Hs=None, CIP=True):
        super().__init__(mol, mol_wo_Hs, CIP)

    def find_fake_axes(self):
        centers = self.find_center_atoms()
        fake_axes = []
        for i in range(len(centers) - 1):
            for j in range(i + 1, len(centers)):
                if self.pub_ring(centers[i], centers[j]) and self.connection[centers[i]][centers[j]]:
                    fake_axes.append((centers[i], centers[j]))
        return fake_axes

    def get_chi_mat(self):
        fake_axes = self.find_fake_axes()
        chi_results = super().get_chi_mat()
        mats, dets, signs = [], [], []
        
        for bond in fake_axes:
            tmp_mat, tmp_det, tmp_sign = [], [], []
            for i in range(len(chi_results["center id"])):
                if bond[0] == chi_results["center id"][i]:
                    tmp_mat.append(chi_results["quadrupole matrix"][i])
                    tmp_sign.append(chi_results["sign"][i])
                    tmp_det.append(chi_results["determinant"][i])
                    break
            for i in range(len(chi_results["center id"])):
                if bond[1] == chi_results["center id"][i]:
                    tmp_mat.append(chi_results["quadrupole matrix"][i])
                    tmp_sign.append(chi_results["sign"][i])
                    tmp_det.append(chi_results["determinant"][i])
                    break
            mats.append(tmp_mat)
            dets.append(tmp_det)
            signs.append(tmp_sign)
        return {"chiral axes": fake_axes, "quadrupole matrix": mats, "determinant": dets, "sign": signs}
