from chiralfinder import ChiralFinder

if __name__ == '__main__':
    smi_list = ["C[C@H]1CC(=O)[C@]2(CCCC2=O)C1", "CC1=CC=C(SC2=C(C)N(C3=CC=CC=C3C(C)(C)C)C(C)=C2)C=C1"]

    chiral_finder = ChiralFinder(smi_list, "SMILES")
    res_ = chiral_finder.get_axial()
    print(res_[0]["chiral axes"], res_[1]["chiral axes"])
    chiral_finder.draw_res_axial("./img_axial", with_index=True)

    smi_list_center = ["BrC/C(=C\[C@@H]1CCCO1)C1CCCCC1"]
    chiral_finder = ChiralFinder(smi_list_center, "SMILES")
    res_ = chiral_finder.get_central()
    print(res_)
