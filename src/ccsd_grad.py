from pyscf.grad import ccsd
from pyscf.lib import logger
from .utilis import relaxed_dm, finalize
from .grad import grad_efield
from .finite_diff import ccsd_grad_fd
import numpy as np
import time

class EFieldCCSDGradients(ccsd.Gradients):
    _keys = ccsd.Gradients._keys
    _keys.update({'efield_strength', 'efield_R', 'efield_rvec', 'efield_atoms','finite_diff'})
    def __init__(self, method):
        super().__init__(method)
        self.efield_strength = method._scf.efield_strength
        self.efield_R = method._scf.efield_R
        self.efield_rvec = method._scf.efield_rvec
        self.efield_atoms = method._scf.efield_atoms
        self.finite_diff = False

    def kernel(self, t1=None, t2=None, l1=None, l2=None, eris=None,
                     atmlst=None, verbose=None):
        log = logger.new_logger(self, verbose)
        mycc = self.base
        if t1 is None:
            if mycc.t1 is None:
                mycc.run()
            t1 = mycc.t1
        if t2 is None: t2 = mycc.t2
        if l1 is None: l1 = mycc.l1
        if l2 is None: l2 = mycc.l2
        if eris is None:
            eris = mycc.ao2mo()
        if t1 is None or t2 is None:
            t1, t2 = mycc.kernel(eris=eris)[1:]
        if l1 is None or l2 is None:
            l1, l2 = mycc.solve_lambda(eris=eris)
        if atmlst is None:
            atmlst = self.atmlst
        else:
            self.atmlst = atmlst

        de = self.grad_elec(t1, t2, l1, l2, eris, atmlst, verbose=log)
        self.de = de + self.grad_nuc(atmlst=atmlst)
        
        # Add electric field contribution
        t0 = time.time()
        dm = relaxed_dm(mycc, self.mol, t1, t2, l1, l2, eris, atmlst)    
        g = grad_efield(self.mol, dm, self.efield_strength, self.efield_rvec, self.efield_R, 
                        mycc._scf.old_paxes, self.efield_atoms, self.finite_diff)
        mycc._scf._set_old_paxes()
        if self.verbose >= logger.NOTE or self.finite_diff:
            finalize(self, g)
            logger.note(self, 'Electric field gradients time: %.3f s' % (time.time()-t0))
        self.de += g

        if self.mol.symmetry:
            self.de = self.symmetrize(self.de, atmlst)
        
        if self.finite_diff:
            de_fd = ccsd_grad_fd(self.mol, self.efield_strength, self.efield_R, self.efield_rvec, self.efield_atoms)
            logger.note(self, 'Check with finite difference %.10f' % np.linalg.norm(de_fd-self.de))
            self._write(self.mol, de_fd, self.atmlst)
        
        self._finalize()
        return self.de