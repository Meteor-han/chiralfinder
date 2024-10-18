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


class ChiralFinder:
    def __init__(self, input_, input_type="SMILES") -> None:
        self.my_classes = [ChiralAxialType1, ChiralAxialType2, ChiralAxialType3, ChiralAxialType4, ChiralAxialType5, ChiralAxialType6]
        self.res_axial = {"v1": [], "v2": [], "v3": [], "v4": [], "v5": [], "v6": [], "merge": []}
        
        self.center_class = ChiralCenter
        self.res_central = []

        self.mols = []

        if input_type == "SMILES":
            for s in input_:
                mol = Chem.MolFromSmiles(s)
                AllChem.EmbedMolecule(mol, maxAttempts=100)
                mol = Chem.AddHs(mol, addCoords=True)
                self.mols.append(mol)
        elif input_type == "molecules":
            for mol in input_:
                if mol.GetNumConformers() == 0:
                    AllChem.EmbedMolecule(mol, maxAttempts=100)
                mol = Chem.AddHs(mol, addCoords=True)
                self.mols.append(mol)
        elif input_type == "sdf":
            r = Chem.SDMolSupplier(input_)
            for mol in tqdm(r):
                if mol.GetNumConformers() == 0:
                    AllChem.EmbedMolecule(mol, maxAttempts=100)
                mol = Chem.AddHs(mol, addCoords=True)
                self.mols.append(mol)
        else:
            warnings.warn("Invalid input format!")
    
    def get_central(self):
        if self.res_central:
            return self.res_central
        for i, mol_original in enumerate(tqdm(self.mols)):
            self.res_central.append(self.center_class(mol_original).get_chi_mat())
        return self.res_central

    def get_axial(self):
        if self.res_axial["merge"]:
            return self.res_axial["merge"]
        for i, mol_original in enumerate(tqdm(self.mols)):
            temp_ = {"chiral axes": [], "quadrupole matrix": [], "determinant": [], "sign": []}

            """merge all res"""
            for j in range(1, 7):
                tag_ = f"v{j}"
                res_o = self.my_classes[j-1](Chem.Mol(mol_original)).get_chi_mat()
                self.res_axial[tag_].append(res_o)
                for i in range(len(res_o["chiral axes"])):
                    if res_o["chiral axes"][i] not in temp_["chiral axes"] and res_o["chiral axes"][i][::-1] not in temp_["chiral axes"]:
                        temp_["chiral axes"].append(res_o["chiral axes"][i])
                        temp_["quadrupole matrix"].append(res_o["quadrupole matrix"][i])
                        temp_["determinant"].append(res_o["determinant"][i])
                        temp_["sign"].append(res_o["sign"][i])
            self.res_axial["merge"].append(temp_)
            
        return self.res_axial["merge"]

    def draw_res_axial(self, dir_path="./img_axial", size=(500, 500)):
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        for i, mol in enumerate(tqdm(self.mols)):
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

            # d = rdMolDraw2D.MolDraw2DSVG(500, 500)
            d = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
            rdMolDraw2D.PrepareAndDrawMolecule(d, mol, highlightAtoms=hit_ats,
                                                highlightAtomColors=atom_cols,
                                                highlightBonds=hit_bonds,
                                                highlightBondColors=bond_cols)
            d.FinishDrawing()

            d.WriteDrawingText(os.path.join(dir_path, f"{i}.png"))

    def draw_res_center(self, dir_path="./img_center", size=(500, 500)):
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        for i, mol in enumerate(tqdm(self.mols)):
            hit_ats = self.res_central[i]["chiral id"]
            colours = [(0., 1.0, 0.)]
            atom_cols = {}
            for j, at in enumerate(hit_ats):
                atom_cols[at] = colours[0]

            # d = rdMolDraw2D.MolDraw2DSVG(500, 500)
            d = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
            rdMolDraw2D.PrepareAndDrawMolecule(d, mol, highlightAtoms=hit_ats,
                                               highlightAtomColors=atom_cols,
                                               )
            d.FinishDrawing()

            d.WriteDrawingText(os.path.join(dir_path, f"{i}.png"))
