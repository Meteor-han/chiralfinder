# chiralfinder

Data and codes for the paper "A Unifying Geometric Framework for Computational Representation of Stereoisomers Based on Mixed Product", in press at *Cell Reports Physical Science*. See pre-print version at [Chemrxiv](https://chemrxiv.org/engage/chemrxiv/article-details/67a0ae146dde43c90873b1e0).

## Overview

Molecular chirality plays a fundamental role in chemistry, biology, and drug discovery, yet existing computational methods struggle with complex stereochemistry and cannot quantify chirality. This project introduces **a unifying geometric framework based on mixed product representations** that directly map molecular symmetry breaking into a 3D algebraic space. The framework provides quantitative descriptors for each stereogenic element, going beyond the qualitative Cahn–Ingold–Prelog rules. We implement this theory in **ChiralFinder**, a computational tool that accurately identifies and distinguishes central and axial stereogenic elements, differentiates conformations, and integrates effectively with machine learning models to improve spectra prediction. By offering a rigorous geometric foundation, ChiralFinder enables automated chirality analysis and forms the basis for a comprehensive stereoisomer representation framework with potential extensions to planar and helical chirality.

<img src="https://github.com/Meteor-han/chiralfinder/blob/main/img_axial/chirality product.png" alt="1"  width="30%" height="auto" />

## Quick use

See online web server at [Here](https://compbio.sjtu.edu.cn/services/chiralfinder).

To use ChiralFinder as a python package, install Anaconda, create and enter your own environment like

    conda create -n env_test python=3.10
Enter the conda environment and install the ChiralFinder package through pip like

```
conda activate env_test
pip install chiralfinder
```

Run `run_example.py` to get example results.

```
python run_example.py
```

```python
from chiralfinder import ChiralFinder

if __name__ == '__main__':
    smi_list = ["C[C@H]1CC(=O)[C@]2(CCCC2=O)C1", "CC1=CC=C(SC2=C(C)N(C3=CC=CC=C3C(C)(C)C)C(C)=C2)C=C1"]

    chiral_finder = ChiralFinder(smi_list, "SMILES")
    res_ = chiral_finder.get_axial(n_cpus=8)
    print(res_[0]["chiral axes"], res_[1]["chiral axes"])
    chiral_finder.draw_res_axial("./img")

    smi_list_center = ["BrC/C(=C\[C@@H]1CCCO1)C1CCCCC1"]
    chiral_finder = ChiralFinder(smi_list_center, "SMILES")
    res_ = chiral_finder.get_central()
    print(res_)
```

You will get the images of two molecules with predicted chiral axes in the folder `./img` by default. Predicted chiral axes:

```
[(5,)] [(9, 10)]
```

<img src="https://github.com/Meteor-han/chiralfinder/blob/main/img_axial/0.png" alt="0" width="30%" height="auto" /><img src="https://github.com/Meteor-han/chiralfinder/blob/main/img_axial/1.png" alt="1"  width="30%" height="auto" />

You will get the prediction of one molecule for central chirality.

```
[{
'center id': [4], 
'quadrupole matrix': 
       [[array([[-0.29989323, -1.08474687,  0.09943544],
       [-2.0754821 ,  0.47857598,  1.02051223],
       [-0.0064714 , -0.03258116,  2.29906673]])]], 
'determinant': [[-5.501797575969392]], 
'norm CP': [[-0.8993660781912431]], 
'sign': [[-1.0]]
}]
```


## Dataset

The RotA dataset is stored in the folder `./data`. The excel file contains labeled chiral axes and some calculated molecular properties. The pickle file includes calculated molecular conformers.

We also provide sampled achiral molecules and centrally chiral molecules with multiple centers from the PubChem3D database in the folder `./data`.

## Citation

```
@article{Shi2025ChiralFinder,
  author    = {Shi, Runhan and Zhang, Chi and Yu, Gufeng and Huo, Xiaohong and Yang, Yang},
  title     = {ChiralFinder: Automated Detection of Stereogenic Elements and Discrimination of Stereoisomers in Complex Molecules},
  journal   = {ChemRxiv},
  year      = {2025},
  doi       = {10.26434/chemrxiv-2025-wz7kh},
  note      = {Preprint, not peer-reviewed}
}
```

