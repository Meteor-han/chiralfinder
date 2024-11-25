from .quadrupole_utils import *
from .quadrupole_center import ChiralCenter
from .quadrupole_v1 import ChiralAxialType1
from .quadrupole_v2 import ChiralAxialType2
from .quadrupole_v3 import ChiralAxialType3
from .quadrupole_v4 import ChiralAxialType4
from .quadrupole_v5 import ChiralAxialType5
from .quadrupole_v6 import ChiralAxialType6
from tqdm import tqdm
from rdkit.Chem.Draw import rdMolDraw2D
import os
from multiprocessing import Pool


axial_classes = [ChiralAxialType1, ChiralAxialType2, ChiralAxialType3, ChiralAxialType4, ChiralAxialType5, ChiralAxialType6]
central_class = ChiralCenter


class ChiralFinder:
    def __init__(self, input_, input_type="SMILES") -> None:
        self.res_axial = {"v1": [], "v2": [], "v3": [], "v4": [], "v5": [], "v6": [], "merge": []}
        self.res_central = []

        self.mols = []

        if input_type == "SMILES":
            for i, s in enumerate(input_):
                mol = Chem.MolFromSmiles(s)
                if not mol:
                    self.mols.append(None)
                    warnings.warn(f"Index: {i}, invalid SMILES: {s}")
                    continue
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, maxAttempts=100)
                self.mols.append(mol)
        elif input_type == "molecules":
            for mol in input_:
                mol = Chem.AddHs(mol)
                if mol.GetNumConformers() == 0:
                    AllChem.EmbedMolecule(mol, maxAttempts=100)
                self.mols.append(mol)
        elif input_type == "sdf":
            for sdf in input_:
                r = Chem.SDMolSupplier(sdf)
                for mol in tqdm(r):
                    mol = Chem.AddHs(mol)
                    if mol.GetNumConformers() == 0:
                        AllChem.EmbedMolecule(mol, maxAttempts=100)
                    self.mols.append(mol)
        elif input_type == "mol":
            for mol_p in input_:
                mol = Chem.MolFromMolFile(mol_p, removeHs=False)
                mol = Chem.AddHs(mol)
                if mol.GetNumConformers() == 0:
                    AllChem.EmbedMolecule(mol, maxAttempts=100)
                self.mols.append(mol)
        else:
            warnings.warn("Invalid input format!")
    
    def get_central(self, n_cpus=4):
        if self.res_central:
            return self.res_central
        with Pool(processes=n_cpus) as pool:
            results = list(tqdm(pool.imap_unordered(self._process_one_mol_central, [(i, mol) for i, mol in enumerate(self.mols)]), total=len(self.mols)))
        pool.close()
        pool.join()

        ordered_results = [result for result in sorted(results, key=lambda x: x[0])]
        for one in ordered_results:
            self.res_central.append(one[1])
        return self.res_central

    def _process_one_mol_central(self, iter_mol):
        index_, mol = iter_mol
        if not mol:
            return index_, None
        return index_, central_class(mol).get_chi_mat()

    def draw_res_center(self, dir_path="./img_center", size=(500, 500), with_index=False):
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        for i, mol in enumerate(tqdm(self.mols)):
            if not mol:
                continue
            # without Hs, remove conformation
            mol = Chem.RemoveHs(mol)
            mol.RemoveAllConformers()
            if with_index:
                mol = mol_with_atom_index(mol)
            hit_ats = self.res_central[i]["center id"]
            colours = [(0., 1.0, 0.)]
            atom_cols = {}
            for j, at in enumerate(hit_ats):
                atom_cols[at] = colours[0]

            d = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
            rdMolDraw2D.PrepareAndDrawMolecule(d, mol, highlightAtoms=hit_ats,
                                               highlightAtomColors=atom_cols,
                                               )
            d.FinishDrawing()

            d.WriteDrawingText(os.path.join(dir_path, f"{i}.png"))
    
    def get_axial(self, n_cpus=4):
        if self.res_axial["merge"]:
            return self.res_axial["merge"]
        with Pool(processes=n_cpus) as pool:
            results = list(tqdm(pool.imap_unordered(self._process_one_mol_axial, [(i, mol) for i, mol in enumerate(self.mols)]), total=len(self.mols)))
        pool.close()
        pool.join()

        ordered_results = [result for result in sorted(results, key=lambda x: x[0])]
        for one in ordered_results:
            self.res_axial["merge"].append(one[1])
            for k, v in one[2].items():
                self.res_axial[k].append(v)
        return self.res_axial["merge"]

    def _process_one_mol_axial(self, iter_mol):
        index_, mol = iter_mol
        if not mol:
            tag2res = {}
            for j in range(1, 7):
                tag2res[f"v{j}"] = None
            return index_, None, tag2res
        temp_ = {"chiral axes": [], "quadrupole matrix": [], "determinant": [], "sign": []}
        tag2res = {}

        """merge all res"""
        for j in range(1, 7):
            tag_ = f"v{j}"
            res_o = axial_classes[j-1](Chem.Mol(mol)).get_chi_mat()
            tag2res[tag_] = res_o
            # self.res_axial[tag_].append(res_o)
            for i in range(len(res_o["chiral axes"])):
                if res_o["chiral axes"][i] not in temp_["chiral axes"] and res_o["chiral axes"][i][::-1] not in temp_["chiral axes"]:
                    temp_["chiral axes"].append(res_o["chiral axes"][i])
                    """for spiral chain, just repeat the res for end atoms, it does not matter"""
                    if len(res_o["quadrupole matrix"]) <= i:
                        i = 0
                    temp_["quadrupole matrix"].append(res_o["quadrupole matrix"][i])
                    temp_["determinant"].append(res_o["determinant"][i])
                    temp_["sign"].append(res_o["sign"][i])
        return index_, temp_, tag2res
    
    def draw_res_axial(self, dir_path="./img_axial", size=(500, 500), with_index=False):
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        for i, mol in enumerate(tqdm(self.mols)):
            if not mol:
                continue
            # without Hs, remove conformation
            mol = Chem.RemoveHs(mol)
            mol.RemoveAllConformers()
            if with_index:
                mol = mol_with_atom_index(mol)
            labels = self.res_axial["merge"][i]["chiral axes"]
            atoms = mol.GetAtoms()
            hit_ats = []
            hit_bonds = []
            for one in labels:
                # check every two
                for k in range(len(one)):
                    for j in range(len(one)):
                        if k != j:
                            atom1 = atoms[one[k]]
                            atom2 = atoms[one[j]]
                            connected = False
                            for bond in atom1.GetBonds():
                                if bond.GetOtherAtomIdx(one[k]) == atom2.GetIdx():
                                    connected = True
                                    hit_bonds.append(mol.GetBondBetweenAtoms(one[k], one[j]).GetIdx())
                                    break
                            if not connected:
                                hit_ats.append(one[k])
                # for single atom
                if len(one) == 1:
                    hit_ats.append(one[0])
            colours = [(0., 1.0, 0.)]
            atom_cols = {}
            for j, at in enumerate(hit_ats):
                atom_cols[at] = colours[0]
            bond_cols = {}
            for j, bd in enumerate(hit_bonds):
                bond_cols[bd] = colours[0]

            d = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
            rdMolDraw2D.PrepareAndDrawMolecule(d, mol, highlightAtoms=hit_ats,
                                                highlightAtomColors=atom_cols,
                                                highlightBonds=hit_bonds,
                                                highlightBondColors=bond_cols)
            d.FinishDrawing()

            d.WriteDrawingText(os.path.join(dir_path, f"{i}.png"))

# add index to the graph
def mol_with_atom_index(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol
