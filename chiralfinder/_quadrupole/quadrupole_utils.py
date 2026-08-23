from rdkit import Chem
from rdkit.Chem import AllChem
from copy import deepcopy
import numpy as np
from collections import defaultdict
import warnings


def normalize_mol_atom_order(mol):
    """Renumber atoms so heavy atoms precede hydrogens.

    RDKit ``AddHs()`` usually appends hydrogens, but molecules built from
    SMILES with explicit H (``removeHs=False``) can place them at arbitrary
    indices.  Downstream code assumes heavy-atom indices match between the
    full molecule and ``Chem.RemoveHs(mol)``.
    """
    if mol is None:
        return None, []

    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1]
    hydrogens = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1]
    new_order = heavy + hydrogens
    if new_order == list(range(mol.GetNumAtoms())):
        return mol, list(range(len(heavy)))
    return Chem.RenumberAtoms(mol, new_order), list(range(len(heavy)))


class ChiralBase:
    def __init__(self, mol=None, mol_wo_Hs=None, max_conf_num=5, sign_eps=1e-4, CIP=True):
        mol, self.heavy_atom_indices = normalize_mol_atom_order(mol)
        self.num_heavy_atoms = len(self.heavy_atom_indices)

        # with Hs
        self.mol = mol
        self.mol_wo_Hs = Chem.RemoveHs(mol) if mol_wo_Hs is None else mol_wo_Hs
        self.atoms = mol.GetAtoms()
        self.ssr = Chem.GetSymmSSSR(mol)
        self.connection = Chem.GetAdjacencyMatrix(mol)
        # coordinates, important; up to max_conf_num
        self.coordinates = []
        self.max_conf_num = max_conf_num
        self.sign_eps = sign_eps
        for conf in mol.GetConformers():
            # np array
            self.coordinates.append(conf.GetPositions())
            if len(self.coordinates) >= max_conf_num:
                break
        # self.conformer = mol.GetConformer()
        self.CIP_list = self._build_cip_list(CIP)

    def _build_cip_list(self, CIP=True):
        canonical_ranks = list(Chem.CanonicalRankAtoms(
            self.mol, breakTies=False, includeChirality=True, includeIsotopes=True))

        if not CIP:
            return canonical_ranks

        cip_list = defaultdict(lambda: -1)
        try:
            # AssignStereochemistry may fail anyway, try AssignStereochemistryFrom3D
            Chem.AssignStereochemistryFrom3D(self.mol_wo_Hs)
            for atom in self.mol_wo_Hs.GetAtoms():
                cip_list[atom.GetIdx()] = int(atom.GetProp('_CIPRank'))
        except Exception:
            warnings.warn("Fail to assign CIPs, use CanonicalRankAtoms instead.")
            return canonical_ranks

        # Heavy-atom CIP ranks live on mol_wo_Hs; fill hydrogen ranks from the
        # full molecule so neighbor lookups never alias a heavy-atom rank.
        for idx in range(self.num_heavy_atoms, self.mol.GetNumAtoms()):
            cip_list[idx] = canonical_ranks[idx]
        return cip_list

    def get_cip_rank(self, atom_idx):
        if isinstance(self.CIP_list, list):
            return self.CIP_list[atom_idx]
        return self.CIP_list[atom_idx]

    # find chiral atoms
    # with Hs, may get false center C "?"
    def find_center_atoms(self):
        c0 = Chem.FindMolChiralCenters(self.mol_wo_Hs, useLegacyImplementation=False)
        cen_chi = [i[0] for i in c0]
        return cen_chi
    
    # find spiral atoms
    def find_spiral_atoms(self):
        # Topology is defined on heavy atoms; indices match mol_wo_Hs after
        # normalize_mol_atom_order().
        ri = self.mol_wo_Hs.GetRingInfo()
        atoms_wo_hs = self.mol_wo_Hs.GetAtoms()
        ssr = Chem.GetSymmSSSR(self.mol_wo_Hs)
        spi_pot = []
        for idx in self.heavy_atom_indices:
            if ri.NumAtomRings(idx) == 2:
                spi_pot.append(idx)

        # delete all atoms who share completely same rings with their neighbors
        spi = deepcopy(spi_pot)
        for j in spi_pot:
            j_ring = []
            for ring in ssr:
                if j in ring:
                    j_ring.append(list(ring))
            if len(j_ring) < 2:
                if j in spi:
                    spi.remove(j)
                continue
            j_ring = j_ring[0] + j_ring[1]

            neighbors = [atom.GetIdx() for atom in atoms_wo_hs[j].GetNeighbors()]
            for k in neighbors:
                if j_ring.count(k) == 2:
                    if j in spi:
                        spi.remove(j)
                    break
        return spi

    # whether two atoms share one ring
    def pub_ring(self, atom_1, atom_2):        
        atom_1_ring = []
        # atom_1 rings
        for ring in self.ssr:
            if atom_1 in list(ring):
                atom_1_ring.extend(list(ring))
        # atom_2
        return True if atom_2 in set(atom_1_ring) else False

    # find 'C=C' like double bonds
    def get_double_bond(self):
        bonds = self.mol.GetBonds()
        double_bonds = []
        for bond in bonds:
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                double_bonds.append(bond)
        return double_bonds

    # whether two bonds share one atom
    def pub_atom(self, bond_1_, bond_2_):
        atom_1_id = [bond_1_.GetBeginAtomIdx(), bond_1_.GetEndAtomIdx()]
        atom_2_id = [bond_2_.GetBeginAtomIdx(), bond_2_.GetEndAtomIdx()]
        criterion = set(atom_1_id + atom_2_id)
        return True if len(criterion) == 3 else False

    # compute the absolute chirality
    def criterion(self, mat):
        det_ = np.linalg.det(mat)
        sign_ = np.sign(det_) if abs(det_) > self.sign_eps else 0
        return det_, sign_

    def check_chain_equal(self, bonds_chain):
        bonds2idx = []
        for one_ in bonds_chain:
            temp = []
            for bond in one_:
                temp.append(sorted([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]))
            bonds2idx.append(sorted(temp))
        return sorted(bonds2idx)

# add Hs, retain existing coor, xy axis reverse coor
def transform_coordinate(mol_):
    # if no conformer, try to embed
    if mol_.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol_, maxAttempts=100)
    mol = Chem.AddHs(mol_, addCoords=True)
    # copy
    mol_r = Chem.Mol(mol)
    mol_r.RemoveAllConformers()
    num_atoms = mol.GetNumAtoms()
    conf_num = mol.GetNumConformers()
    assert conf_num > 0.5
    for i in range(conf_num):
        c = mol.GetConformer(i)
        p = c.GetPositions()
        p_r = np.copy(p)
        p_r[:, 2] = -p_r[:, 2]
        _conformer = Chem.Conformer(num_atoms)
        for j in range(len(p_r)):
            _conformer.SetAtomPosition(j, p_r[j])
        mol_r.AddConformer(_conformer, assignId=True)

    return mol, mol_r
