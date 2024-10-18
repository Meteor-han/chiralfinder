from quadrupole_utils import *


class ChiralAxialType5(ChiralBase):
    def __init__(self, mol):
        super().__init__(mol)

    # 寻找不在环内的单键
    def get_single_bonds(self):
        # 获取单键模板
        # A = Chem.MolFromSmiles('CC')
        # bonds_A = A.GetBonds()
        # single = bonds_A[0].GetBondType()
        single = Chem.BondType.SINGLE

        # 获取不在环内的单键
        bonds = self.mol.GetBonds()  # if self.mol, may get -H, too many; but we need Hs right?
        single_bonds = []
        ring_bonds=[]
        for bond in bonds:
            if bond.GetBondType() == single:  # 与单键模板比较
                for i in range(1, 13):  #12元环以下为中小环，单键旋转受阻
                    if(bond.IsInRingSize(i)):  #若单键在中小环里，则保留
                        ring_bonds.append(bond)
                        break

                atom_1, atom_2 = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())  # 获取单键两端原子
                if not self.pub_ring(atom_1, atom_2):  # 判断单键两端原子是否共环，这等价于判断单键是否在某个环内
                    single_bonds.append(bond)
        return single_bonds, ring_bonds  # 返回所有不在环内的单键构成的列表，其元素为键对象， 返回分子中所有在中小环里单键的列表

    '''
    对于某一根不在环内的单键，将分子中的剩余原子分成两类
    一类与atom_1相连，一类与atom_2相连
    '''
    def classify_atoms(self, bond, n):
        assert n > 0
        # mol = Chem.AddHs(mol)
        atoms = self.mol.GetAtoms()  # 获取分子的原子信息
        atom_1, atom_2 = (bond.GetBeginAtom(), bond.GetEndAtom())  # 获取键两端原子
        
        nei_1 = set([atom_1.GetIdx()])
        for _ in range(n):  # 迭代次数预先设定，只考虑以atom_1为“圆心”，n个原子为“半径”的范围（类似摩根指纹的想法）
            nei_temp = set()
            for atom in nei_1:  # 获取邻居
                nei = atoms[atom].GetNeighbors()
                nei_temp.update([atom.GetIdx() for atom in nei])
            nei_1.update(nei_temp)

            # 掐灭往atom_2方向生长的趋势
            nei_1.discard(atom_2.GetIdx())
            # 除去atom_1，避免后面程序中计算余弦时出现问题
            nei_1.discard(atom_1.GetIdx())

        # 获取n个原子以内相连的原子-atom_2（与atom_1基本一致）
        nei_2 = set([atom_2.GetIdx()])
        for _ in range(n):
            nei_temp = set()
            for atom in nei_2:  # 获取邻居
                nei = atoms[atom].GetNeighbors()
                nei_temp.update([atom.GetIdx() for atom in nei])
            nei_2.update(nei_temp)

            # 掐灭往atom_1方向生长的趋势
            nei_2.discard(atom_1.GetIdx())
            # 除去atom_2，避免后面程序中计算余弦时出现问题
            nei_2.discard(atom_2.GetIdx())

        return [nei_1, nei_2]  # 返回一个二维列表

    '''
    计算同侧共面之后的平面距离并以此判断是否阻转
    '''
    def calculate_planar_distance(self, nei_1, nei_2, atom_1, atom_2, conf_id):
        """
        原子半径数据，单位为pm,，乘0.01得以埃为单位的半径，与mol文件中一致
        来源《无机化学》（上册） 北京师范大学 华中师范大学 南京师范大学编
        前五周期的元素的原子半径，其中金属元素采用金属半径，非金属元素采用共价半径
        """
        # atom_r = np.array([30, 140, 152, 111.3, 88, 77.2, 70, 66, 64, 154, 186, 160, 143.1, 117, 110, 104, 99, 192,
        #                 232, 197, 162, 147, 134, 128, 127, 126, 125, 124, 128, 134, 135, 128, 121, 117, 114, 198,
        #                 248, 215, 180, 160, 146, 139, 136, 134, 134, 137, 144, 148.9, 167, 151, 145, 137, 133,
        #                 218]) * 0.01
        atom_r = np.array([30,140,152,111.3,88,77.2,70,66,64,154,186,160,143.1,117,110,104,99,192,
                            232,197,162,147,134,128,127,126,125,124,128,134,135,128,121,117,114,198,
                            248,215,180,160,146,139,136,134,134,137,144,148.9,167,151,145,137,133,218,
                            265,217.3,183,181.8,182.4,183.4,180.4,208.4,180.4,177.3,178.1,176.2,176.1,175.9,193.3,173.8,
                            159,146,139,137,135,135.5,138.5,144,151,170,175,154.7,164])*0.01

        atoms = self.mol.GetAtoms()  # 获取原子信息
        atom_1, atom_2 = atom_1.GetIdx(), atom_2.GetIdx()
        atom_1_cor = self.coordinates[conf_id][atom_1]  # 获取atom_1的坐标信息，作为原点
        atom_2_cor = self.coordinates[conf_id][atom_2] - atom_1_cor  # 计算atom_2在以atom_1为原点的坐标系中的坐标
        eps = 0.001  # 初始设定的一个小数，用于判断判别式是否足够小
        t = 0  # 计数器，用来判定是否阻转
        for i in nei_1:
            for j in nei_2:
                #获取坐标信息
                atom_i_cor = self.coordinates[conf_id][i] - atom_1_cor  # 计算atom_i在以atom_1为原点的坐标系中的坐标
                atom_j_cor = self.coordinates[conf_id][j] - atom_1_cor  # 计算atom_j在以atom_1为原点的坐标系中的坐标
                
                #转化为极坐标
                #到原点的距离
                r_i = np.linalg.norm(atom_i_cor)
                r_j = np.linalg.norm(atom_j_cor)
                r_x = np.linalg.norm(atom_2_cor)
                
                #与极轴的夹角
                cos_i = np.clip(np.dot(atom_i_cor, atom_2_cor) / (r_i*r_x), -1., 1.)
                cos_j = np.clip(np.dot(atom_j_cor, atom_2_cor) / (r_j*r_x), -1., 1.)
                
                #计算旋转后同侧平面距离
                cos_i_j = np.cos(np.arccos(cos_i) - np.arccos(cos_j))
                plan_dis = np.sqrt(r_i**2 + r_j**2 - 2*r_i*r_j*cos_i_j)

                # 根据原子序数获取原子半径
                r_1 = atom_r[atoms[i].GetAtomicNum() - 1]
                r_2 = atom_r[atoms[j].GetAtomicNum() - 1]

                # 计算判别式的值并判断
                delta = plan_dis - r_1 - r_2
                if delta < eps:  # 若判别式值足够小，则为阻转轴
                    t += 1
                    break
        
        return True if t > 0 else False

    # 主程序
    def get_chi_mat(self, n=10):  # mol_0:输入带有坐标信息的分子，n:考虑n个原子以内的邻居
        # mol=Chem.MolFromSmiles(mol_0)#读取分子
        bonds, ring_bonds = self.get_single_bonds()  # 获取所有环外单键, #寻找所有在中小环里的单键
        rotation_limited_idx = [[] for _ in range(len(self.coordinates))]
        rotation_limited = set()
        for bond in bonds:  # 遍历所有环外单键
            atom_1, atom_2 = (bond.GetBeginAtom(), bond.GetEndAtom())  # 获取键的两端原子
            classify = self.classify_atoms(bond, n)  # 把分子中的原子按atom_1,atom_2分类
            # if bond.GetBeginAtomIdx()==3 or bond.GetBeginAtomIdx()==4:
            #     print()
            for id_ in range(len(self.coordinates)):
                if self.calculate_planar_distance(classify[0], classify[1], atom_1, atom_2, id_):
                    rotation_limited_idx[id_].append((atom_1.GetIdx(), atom_2.GetIdx()))
                    rotation_limited.add(bond)
        # return rotation_limited  # 返回所有阻转轴（单键），如果没有阻转轴，则为空列表
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

    #增加部分：寻找环内阻转轴及进一步寻找手性轴
    #在阻转轴里寻找手性轴
    def find_chiral_axes(self, rotation_limited, ring_bonds):
        # mol = Chem.AddHs(mol)#分子加氢;已经有H了
        # atoms = self.mol.GetAtoms()  #分子中所有原子的信息
        # just an order, not CIP, unable to relate to R/S
        res = list(Chem.CanonicalRankAtoms(self.mol, breakTies=False, includeChirality=True, includeIsotopes=True))
        bonds = self.bonds_set(rotation_limited + ring_bonds)  #合并“阻转轴”，其实已经去重了，不需要bonds_set

        #获取手性轴
        chiral_axes = []  # merge all confs
        mats, dets, signs = [], [], []  # for each conf
        for bond in bonds:
            #获取阻转轴两端原子的邻居
            begin_neighbor = [atom.GetIdx() for atom in bond.GetBeginAtom().GetNeighbors()]
            end_neighbor = [atom.GetIdx() for atom in bond.GetEndAtom().GetNeighbors()]
            #删去邻居中所包含的另一端原子，得到界外邻居
            if bond.GetEndAtomIdx() in begin_neighbor:
                begin_neighbor.remove(bond.GetEndAtomIdx())
            if bond.GetBeginAtomIdx() in end_neighbor:
                end_neighbor.remove(bond.GetBeginAtomIdx())
        
            #大部分轴手性每一端界外邻居为2个，仅考虑此种情况
            if (len(begin_neighbor)!=2) or (len(end_neighbor)!=2):
                continue

           #若每一端的2个界外邻居均相异，则具有轴手性
            if (res[begin_neighbor[0]] != res[begin_neighbor[1]]) and (res[end_neighbor[0]] != res[end_neighbor[1]]):  
                chiral_axes.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
                
                mat_confs = []
                det_confs = []
                sign_confs = []
                for conf_ in self.coordinates:
                    begin_cor = conf_[bond.GetBeginAtomIdx()]
                    end_cor = conf_[bond.GetEndAtomIdx()]
                
                    #对轴端邻居进行排序
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
                    
                    #计算手性矩阵
                    # print(begin_neighbor, end_neighbor)
                    a = begin_1_cor - (begin_cor+end_cor)/2
                    b = begin_2_cor - (begin_cor+end_cor)/2
                    c = end_1_cor - end_2_cor
                    mat = np.array([a, b, c])
                    
                    # if (abs(np.linalg.det(A))>eps):
                    mat_confs.append(mat)
                    det_, sign_ = self.criterion(mat)
                    det_confs.append(det_)
                    sign_confs.append(sign_)
                mats.append(mat_confs)
                dets.append(det_confs)
                signs.append(sign_confs)
        return {"chiral axes": chiral_axes, "quadrupole matrix": mats, "determinant": dets, "sign": signs}
        # return chiral_axes, chiral_mat  #返回阻转轴
