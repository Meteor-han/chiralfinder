from .quadrupole_utils import *


class ChiralAxialType5(ChiralBase):
    def __init__(self, mol):
        super().__init__(mol)

    def get_single_bonds(self):
        # template
        single = Chem.BondType.SINGLE
    
        bonds = self.mol.GetBonds()
        single_bonds = []
        ring_bonds=[]
        for bond in bonds:
            if bond.GetBondType() == single:
                for i in range(1, 13):
                    if(bond.IsInRingSize(i)):
                        ring_bonds.append(bond)
                        break

                atom_1, atom_2 = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                if not self.pub_ring(atom_1, atom_2): 
                    single_bonds.append(bond)
        return single_bonds, ring_bonds 

    def classify_atoms(self, bond, n):
        assert n > 0
        atoms = self.mol.GetAtoms() 
        atom_1, atom_2 = (bond.GetBeginAtom(), bond.GetEndAtom()) 
        
        nei_1 = set([atom_1.GetIdx()])
        for _ in range(n): 
            nei_temp = set()
            for atom in nei_1:
                nei = atoms[atom].GetNeighbors()
                nei_temp.update([atom.GetIdx() for atom in nei])
            nei_1.update(nei_temp)

            nei_1.discard(atom_2.GetIdx())
            nei_1.discard(atom_1.GetIdx())

        nei_2 = set([atom_2.GetIdx()])
        for _ in range(n):
            nei_temp = set()
            for atom in nei_2: 
                nei = atoms[atom].GetNeighbors()
                nei_temp.update([atom.GetIdx() for atom in nei])
            nei_2.update(nei_temp)

            nei_2.discard(atom_1.GetIdx())
            nei_2.discard(atom_2.GetIdx())

        return [nei_1, nei_2]

    def calculate_planar_distance(self, nei_1, nei_2, atom_1, atom_2, conf_id):
        atom_r = np.array([30,140,152,111.3,88,77.2,70,66,64,154,186,160,143.1,117,110,104,99,192,
                            232,197,162,147,134,128,127,126,125,124,128,134,135,128,121,117,114,198,
                            248,215,180,160,146,139,136,134,134,137,144,148.9,167,151,145,137,133,218,
                            265,217.3,183,181.8,182.4,183.4,180.4,208.4,180.4,177.3,178.1,176.2,176.1,175.9,193.3,173.8,
                            159,146,139,137,135,135.5,138.5,144,151,170,175,154.7,164])*0.01

        atoms = self.mol.GetAtoms() 
        atom_1, atom_2 = atom_1.GetIdx(), atom_2.GetIdx()
        atom_1_cor = self.coordinates[conf_id][atom_1] 
        atom_2_cor = self.coordinates[conf_id][atom_2] - atom_1_cor 
        eps = 0.001 
        t = 0 
        for i in nei_1:
            for j in nei_2:
                atom_i_cor = self.coordinates[conf_id][i] - atom_1_cor
                atom_j_cor = self.coordinates[conf_id][j] - atom_1_cor
                
                r_i = np.linalg.norm(atom_i_cor)
                r_j = np.linalg.norm(atom_j_cor)
                r_x = np.linalg.norm(atom_2_cor)
                
                cos_i = np.clip(np.dot(atom_i_cor, atom_2_cor) / (r_i*r_x), -1., 1.)
                cos_j = np.clip(np.dot(atom_j_cor, atom_2_cor) / (r_j*r_x), -1., 1.)
                
                cos_i_j = np.cos(np.arccos(cos_i) - np.arccos(cos_j))
                plan_dis = np.sqrt(r_i**2 + r_j**2 - 2*r_i*r_j*cos_i_j)

                r_1 = atom_r[atoms[i].GetAtomicNum() - 1]
                r_2 = atom_r[atoms[j].GetAtomicNum() - 1]

                delta = plan_dis - r_1 - r_2
                if delta < eps:
                    t += 1
                    break
        
        return True if t > 0 else False

    def get_chi_mat(self, n=10):
        bonds, ring_bonds = self.get_single_bonds() 
        rotation_limited_idx = [[] for _ in range(len(self.coordinates))]
        rotation_limited = set()
        for bond in bonds: 
            atom_1, atom_2 = (bond.GetBeginAtom(), bond.GetEndAtom()) 
            classify = self.classify_atoms(bond, n) 
            for id_ in range(len(self.coordinates)):
                if self.calculate_planar_distance(classify[0], classify[1], atom_1, atom_2, id_):
                    rotation_limited_idx[id_].append((atom_1.GetIdx(), atom_2.GetIdx()))
                    rotation_limited.add(bond)
        return self.find_chiral_axes(list(rotation_limited), ring_bonds)

    def bonds_set(self, bonds):
        set_bonds_idx = set()
        list_bonds = []
        for bond in bonds:
            bond_idx = (bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx())
            if bond_idx not in set_bonds_idx:
                set_bonds_idx.add(bond_idx)
                list_bonds.append(bond)
        return list_bonds

    def find_chiral_axes(self, rotation_limited, ring_bonds):
        res = list(Chem.CanonicalRankAtoms(self.mol, breakTies=False, includeChirality=True, includeIsotopes=True))
        bonds = self.bonds_set(rotation_limited + ring_bonds)

        chiral_axes = []  # merge all confs
        mats, dets, signs = [], [], []  # for each conf
        for bond in bonds:
            begin_neighbor = [atom.GetIdx() for atom in bond.GetBeginAtom().GetNeighbors()]
            end_neighbor = [atom.GetIdx() for atom in bond.GetEndAtom().GetNeighbors()]
            if bond.GetEndAtomIdx() in begin_neighbor:
                begin_neighbor.remove(bond.GetEndAtomIdx())
            if bond.GetBeginAtomIdx() in end_neighbor:
                end_neighbor.remove(bond.GetBeginAtomIdx())
        
            if (len(begin_neighbor)!=2) or (len(end_neighbor)!=2):
                continue

            if (res[begin_neighbor[0]] != res[begin_neighbor[1]]) and (res[end_neighbor[0]] != res[end_neighbor[1]]):  
                chiral_axes.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
                
                mat_confs = []
                det_confs = []
                sign_confs = []
                for conf_ in self.coordinates:
                    begin_cor = conf_[bond.GetBeginAtomIdx()]
                    end_cor = conf_[bond.GetEndAtomIdx()]
                
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
                    
                    a = begin_1_cor - (begin_cor+end_cor)/2
                    b = begin_2_cor - (begin_cor+end_cor)/2
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
