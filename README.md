# chiralfinder

Data and codes for the paper [To be filled].

## Quick use

Install Anaconda, create and enter your own environment like

    conda create -n env_test python=3.10
Enter the conda environment and install the ChiralFinder package through pip like

```
conda activate env_test
pip install chiralfinder
```

Run `test.py` to get example results.

```
python test.py
```

```python
from chiralfinder import ChiralFinder

if __name__ == '__main__':
    smi_list = ["C[C@H]1CC(=O)[C@]2(CCCC2=O)C1", "CC1=CC=C(SC2=C(C)N(C3=CC=CC=C3C(C)(C)C)C(C)=C2)C=C1"]

    chiral_finder = ChiralFinder(smi_list, "SMILES")
    res_ = chiral_finder.get_axial()
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

<img src="./img_axial/0.png" alt="0" style="zoom:40%;" /><img src="./img_axial/1.png" alt="1" style="zoom:40%;" />

You will get the prediction of one molecule for central chirality.

```
[{
'center id': [4], 
'quadrupole matrix': 
       [[array([[-0.29989323, -1.08474687,  0.09943544],
       [-2.0754821 ,  0.47857598,  1.02051223],
       [-0.0064714 , -0.03258116,  2.29906673]])]], 
'determinant': [[-5.501797575969392]], 
'sign': [[-1.0]]
}]
```




## The RotA dataset

The RotA dataset in excel with labeled chiral axes and some calculated molecular properties is stored in the folder `./data`.

## Citation

```
To be filled
```

