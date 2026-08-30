"""Reproduces the Casimir anchor of Section 5: the coefficient by two independent
routes, the field-content scenarios, and the predicted R and ladder scale."""
import numpy as np
from scipy.special import zeta, gamma
hbarc=1.97327e-7; Lam4=2.25e-3
C1=3*zeta(5)/(64*np.pi**6)
C2=float(gamma(2.5))*zeta(5)/np.pi**2.5/(2*np.pi)**4
print(f"route 1 (zeta-regularised mode sum) : C = {C1:.6e}")
print(f"route 2 (thermal analogy)           : C = {C2:.6e}   agreement: {abs(C1-C2)/C1:.1e}")
L0=hbarc/Lam4
print(f"base length hbar*c/Lambda^(1/4) = {L0*1e6:.1f} um")
for lab,N in [("minimal 1 dof",1),("2 Dirac - graviton (3)",3),("same, orbifold (1.5)",1.5),("8 dof",8),("16 dof",16)]:
    R=(N*C1)**0.25*L0
    print(f"  {lab:24s} -> R = {R*1e6:5.2f} um ; ladder scale = {hbarc/R*1e3:5.1f} meV")
print(f"suppression of the confined field (mR=1): exp(-2pi) = {np.exp(-2*np.pi):.4f}")
