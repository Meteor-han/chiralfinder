import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from chiralfinder._quadrupole.quadrupole_utils import normalize_mol_atom_order
from chiralfinder._quadrupole.quadrupole_v1 import ChiralAxialType1
from chiralfinder._quadrupole.quadrupole_center import ChiralCenter


def _embed(mol, seed=7):
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    return mol


def _mol_with_h_at_front_and_middle(base_mol):
    heavy = [a.GetIdx() for a in base_mol.GetAtoms() if a.GetAtomicNum() != 1]
    hydrogens = [a.GetIdx() for a in base_mol.GetAtoms() if a.GetAtomicNum() == 1]
    if len(hydrogens) < 2:
        pytest.skip("molecule does not have enough hydrogens for this permutation")

    mid = len(heavy) // 2
    new_order = [hydrogens[0]] + heavy[:mid] + [hydrogens[1]] + heavy[mid:] + hydrogens[2:]
    return Chem.RenumberAtoms(base_mol, new_order)


def test_issue_2_reproduction_h_at_index_zero():
    """Regression for https://github.com/Meteor-han/chiralfinder/issues/2"""
    smi = "[H]c1ccc2ccccc2c1"
    params = Chem.SmilesParserParams()
    params.removeHs = False
    mol = Chem.MolFromSmiles(smi, params)
    Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL, catchErrors=True)
    mol = _embed(Chem.AddHs(mol))

    assert mol.GetAtomWithIdx(0).GetSymbol() == "H"

    res = ChiralAxialType1(Chem.Mol(mol)).get_chi_mat()
    assert "chiral axes" in res
    assert "sign" in res


def test_hydrogen_in_the_middle_of_atom_list():
    smi = "[H]c1ccc2ccccc2c1"
    params = Chem.SmilesParserParams()
    params.removeHs = False
    mol = Chem.MolFromSmiles(smi, params)
    Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL, catchErrors=True)
    mol = _embed(Chem.AddHs(mol))

    heavy_count = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() != 1)
    mol = _mol_with_h_at_front_and_middle(mol)
    h_in_middle = 1 + heavy_count // 2

    assert mol.GetAtomWithIdx(0).GetSymbol() == "H"
    assert mol.GetAtomWithIdx(h_in_middle).GetSymbol() == "H"

    res = ChiralAxialType1(Chem.Mol(mol)).get_chi_mat()
    assert "chiral axes" in res


def test_normalize_mol_atom_order_aligns_heavy_indices():
    smi = "[H]c1ccc2ccccc2c1"
    params = Chem.SmilesParserParams()
    params.removeHs = False
    mol = Chem.MolFromSmiles(smi, params)
    mol = _embed(Chem.AddHs(mol))

    normalized, heavy_indices = normalize_mol_atom_order(mol)
    mol_wo_hs = Chem.RemoveHs(normalized)

    for idx in heavy_indices:
        assert normalized.GetAtomWithIdx(idx).GetAtomicNum() != 1
        assert mol_wo_hs.GetAtomWithIdx(idx).GetAtomicNum() != 1

    for idx in range(len(heavy_indices), normalized.GetNumAtoms()):
        assert normalized.GetAtomWithIdx(idx).GetAtomicNum() == 1


def test_central_chirality_with_h_at_index_zero():
    smi = "C[C@H](O)Cl"
    params = Chem.SmilesParserParams()
    params.removeHs = False
    mol = Chem.MolFromSmiles(smi, params)
    Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL, catchErrors=True)
    mol = _embed(Chem.AddHs(mol))

    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1]
    hs = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1]
    mol = Chem.RenumberAtoms(mol, hs[:1] + heavy + hs[1:])

    assert mol.GetAtomWithIdx(0).GetSymbol() == "H"

    res = ChiralCenter(Chem.Mol(mol)).get_chi_mat()
    assert "center id" in res
