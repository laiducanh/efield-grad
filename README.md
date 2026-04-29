Analytic nuclear gradients in the presence of oriented external electric fields in a molecule-fixed frame.

This repository provides implementations for single-point calculations and geometry optimizations under oriented external electric fields within molecule-fixed frames, including Principal axis frame and Local reference frame.

---
### Requirements
- Python >= 3.9
- PySCF
- NumPy
- SciPy

### Installation
1. Clone the repository
```bash
git clone https://github.com/laiducanh/efield-grad.git
cd efield-grad
```
2. Install dependencies
It is recommended to use a virtual environment:
```bash
python -m venv efield
source efield/bin/activate
```
Install PySCF (required):
```bash
pip install pyscf
```
Install additional dependencies:
```bash
pip install numpy scipy
```

### Usage Example
First, build a PySCF molecule instance
```Python
import numpy as np
from pyscf import gto

mol = gto.Mole()
mol.atom = '''
  H    1.4353757   -0.8616269   -0.4918772
  C    0.3807982   -1.0704377   -0.2704816
  H   -0.2230878   -0.5142256   -0.9996510
  H    0.2090717   -2.1404517   -0.4428810
  C    0.0239969   -0.6961191    1.1639050
  H    0.6592733   -1.2830607    1.8725512
  O    0.3398130    0.6769092    1.3173146
  H    0.0422517    0.9291206    2.1810295
  C   -1.4501239   -0.9786682    1.4492349
  H   -2.1113592   -0.3953701    0.7946739
  H   -1.6792480   -2.0395972    1.2886605
  H   -1.7184479   -0.7378270    2.4864575
'''
mol.basis = 'cc-pvdz'
mol.unit = 'Angstrom'
mol.build()
```

