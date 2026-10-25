"""Regression coverage for issue #3 and small batches of repository data."""
import pickle
from pathlib import Path
import zipfile

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdCIPLabeler

from chiralfinder import ChiralFinder
from chiralfinder._quadrupole.run import axial_classes
from chiralfinder._quadrupole.quadrupole_center import ChiralCenter
from chiralfinder._quadrupole.quadrupole_utils import ChiralBase, normalize_mol_atom_order


DATA = Path(__file__).resolve().parents[1] / "data"


def _interleave_hydrogens(mol):
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1]
    hs = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1]
    order = []
    for i in range(max(len(heavy), len(hs))):
        order.extend(hs[i:i + 1])
        order.extend(heavy[i:i + 1])
    return Chem.RenumberAtoms(mol, order)


def _mirror(mol):
    mol = Chem.Mol(mol)
    for conf in mol.GetConformers():
        for i, pos in enumerate(conf.GetPositions()):
            pos[2] *= -1
            conf.SetAtomPosition(i, pos)
    return mol


def _assert_numeric_equal(left, right, factor=1):
    if isinstance(left, (list, tuple)):
        assert len(left) == len(right)
        for a, b in zip(left, right):
            _assert_numeric_equal(a, b, factor)
    else:
        np.testing.assert_allclose(left, np.asarray(right) * factor, atol=1e-10)


def _check_geometry(cls, mol, identity):
    original = cls(Chem.Mol(mol))
    result = original.get_chi_mat()
    shuffled = cls(_interleave_hydrogens(mol)).get_chi_mat()
    reflected = cls(_mirror(mol))
    # Hold priority order fixed to test geometry: pseudoasymmetric ligands can
    # exchange CIP priorities when the entire molecule is reflected.
    reflected.CIP_list = original.CIP_list
    mirrored = reflected.get_chi_mat()
    # Some axial branches collect axes through sets; output order is not stable.
    for output in (result, shuffled, mirrored):
        order = sorted(range(len(output[identity])), key=lambda i: output[identity][i])
        for key in (identity, "determinant", "sign"):
            # v2 shares one descriptor across the atoms of a spiro chain.
            if len(output[key]) == len(order):
                output[key] = [output[key][i] for i in order]
    assert result[identity] == shuffled[identity] == mirrored[identity]
    for key in ("determinant", "sign"):
        _assert_numeric_equal(result[key], shuffled[key])
        _assert_numeric_equal(result[key], mirrored[key], factor=-1)
    return result


def test_issue_3_public_api():
    smiles = ["C[C@H](O)CC", "C[C@@H](O)CC", "F[C@H](Cl)Br",
              "C[C@](O)(CC)c1ccccc1"]
    mols = [Chem.AddHs(Chem.MolFromSmiles(s)) for s in smiles]
    for mol in mols:
        assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    results = ChiralFinder(mols, input_type="molecules").get_central(n_cpus=1)
    for mol, result, sign in zip(mols, results, [1, -1, -1, 1]):
        assert result["sign"] == [[sign]]
        neighbors = result["neighbor ids"][0]
        if any(mol.GetAtomWithIdx(i).GetAtomicNum() == 1 for i in neighbors):
            assert mol.GetAtomWithIdx(neighbors[0]).GetAtomicNum() == 1
        _check_geometry(ChiralCenter, mol, "center id")


def test_canonical_mode_and_fallback_unchanged(monkeypatch):
    mol = Chem.AddHs(Chem.MolFromSmiles("C[C@H](O)CC"))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    expected = list(Chem.CanonicalRankAtoms(mol, breakTies=False,
                    includeChirality=True, includeIsotopes=True))
    assert ChiralBase(Chem.Mol(mol), CIP=False).CIP_list == expected

    def fail(*args, **kwargs):
        raise ValueError("forced assignment failure")

    monkeypatch.setattr(Chem, "AssignStereochemistryFrom3D", fail)
    with pytest.warns(UserWarning, match="Fail to assign CIPs"):
        assert ChiralBase(Chem.Mol(mol)).CIP_list == expected


def test_central_dataset_small_batch():
    with zipfile.ZipFile(DATA / "Multi_center.zip") as archive:
        with archive.open("Multi_center.pkl") as stream:
            records = list(pickle.load(stream).values())
    centers = 0
    non_rs = 0
    mismatches = []
    for idx in np.linspace(0, len(records) - 1, 24, dtype=int):
        mol = Chem.AddHs(records[idx][0][0], addCoords=True)
        mol, _ = normalize_mol_atom_order(mol)
        result = _check_geometry(ChiralCenter, mol, "center id")
        labeled = Chem.RemoveHs(mol)
        Chem.AssignStereochemistryFrom3D(labeled)
        Chem.FindMolChiralCenters(labeled, useLegacyImplementation=False)
        rdCIPLabeler.AssignCIPLabels(labeled)
        assert result["center id"]
        for center, signs in zip(result["center id"], result["sign"]):
            atom = labeled.GetAtomWithIdx(center)
            label = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None
            # Pseudoasymmetric/unassigned centers have no simple R/S sign oracle;
            # their permutation and mirror behavior is still checked above.
            if label not in ("R", "S"):
                non_rs += 1
                continue
            if signs != [{"S": 1, "R": -1}[label]]:
                mismatches.append((int(idx), center, label, signs))
            centers += 1
    assert centers > 0
    print(f"central batch: 24 molecules, {centers} R/S centers, {non_rs} other centers")
    # Verified unchanged against c3c653a^ (before issue #2's ranking change).
    # The legacy _CIPRank order is not a complete implementation of modern CIP.
    assert mismatches == [(14019, 2, "R", [1.0]), (16971, 4, "R", [1.0])]


def test_axial_public_api():
    smiles = ["C[C@H]1CC(=O)[C@]2(CCCC2=O)C1",
              "CC1=CC=C(SC2=C(C)N(C3=CC=CC=C3C(C)(C)C)C(C)=C2)C=C1"]
    mols = [Chem.AddHs(Chem.MolFromSmiles(s)) for s in smiles]
    for mol in mols:
        assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    results = ChiralFinder(mols, input_type="molecules").get_axial(n_cpus=1)
    assert [r["chiral axes"] for r in results] == [[(5,)], [(9, 10)]]


def test_axial_dataset_small_batch():
    with (DATA / "RotA.pkl").open("rb") as stream:
        mols = list(pickle.load(stream).values())
    counts = [0] * 6
    conformers = 0
    for idx in np.linspace(0, len(mols) - 1, 36, dtype=int):
        mol = Chem.AddHs(mols[idx], addCoords=True)
        conformers += mol.GetNumConformers()
        for i, cls in enumerate(axial_classes):
            result = _check_geometry(cls, mol, "chiral axes")
            counts[i] += len(result["chiral axes"])
    assert all(counts), counts
    print(f"axial batch: 36 molecules, {conformers} conformers, axes by v1-v6: {counts}")
