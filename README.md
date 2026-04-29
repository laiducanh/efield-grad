Analytic nuclear gradients in the presence of oriented external electric fields in a molecule-fixed frame.

This repository provides implementations for **single-point calculations** and **geometry optimizations** under oriented external electric fields within molecule-fixed frames, including:

- Principal Axis Frame (PAF)
- Local Reference Frame (LRF)

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
1. Build a PySCF molecule
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
2. Define the electric field

The electric field vector $E$ is defined as 

$$
E = U R r \Vert E\Vert
$$

where $U$, $R$, $r$, and $\Vert E\Vert$ are the transformation matrix, rotation matrix, reference vector, and field strength (in atomic units), respectively. The transformation matrix $U$ will be internally determined depending on the molecule-fixed frame in used. By default, the electric field is assumed to be defined in the PAF. In this frame, the eigenvectors of the inertia tensor are sorted in descending order, such that the first principal axis corresponds to the largest eigenvalue. If three non-collinear atoms are specified, the LRF is used instead.
```Python
from scipy.spatial.transform import Rotation

efield = 0.05
R = Rotation.from_euler('x', 90, degrees=True).as_matrix()
rvec=[0, 0, 1]
# specify three atoms for the LRF
atoms=[5, 2, 6]
# atoms = None # otherwise the PAF is used by default
```
3. Single-point calculations

Restricted Hartree-Fock and CCSD single-point calculations are performed with field-perturbed Hamiltonian
```Python
from efield_grad import EFieldRHF

mf = EFieldRHF(mol, efield=efield, R=R, rvec=rvec, atoms=atoms)
mf.kernel()
mycc = cc.CCSD(mf)
# mycc.set_frozen() # frozen-core
e_corr, t1, t2 = mycc.kernel()
```
4. Geometry optimization

Geometry optimization using CCSD level
```Python
from efield_grad import EFieldCCSDGradients

conv_params = { # These are the default settings
    'convergence_energy': 1e-5,  # Eh
    'convergence_grms': 3e-4,    # Eh/Bohr
    'convergence_gmax': 4.5e-4,  # Eh/Bohr
    'convergence_drms': 1.2e-3,  # Angstrom
    'convergence_dmax': 1.8e-3,  # Angstrom
    'maxiter': 1,    
}
mol_eq = EFieldCCSDGradients(mycc)
mol_eq.optimizer(solver='geometric').kernel(conv_params)
print(mol_eq.tostring())
```

### Citation
If you use this code in your research, please cite:
> Duc Anh Lai, Devin Matthews (2026). *Analytic nuclear gradients including oriented external electric fields in a molecule-fixed frame*. J. Chem. Theory Comput. (under review) [arXiv:2604.01189](https://arxiv.org/abs/2604.01189).