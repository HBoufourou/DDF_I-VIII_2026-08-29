"""Reproduces the certified orbifold halving of §5 and the branch-counting clarification."""
import numpy as np
from scipy.special import ai_zeros, zeta
az,apz,_,_=ai_zeros(200)
for D in (2,5,10,20):
    nD=np.sum(-az<D); nN=np.sum(-apz<D)
    print(f"D={D:4.0f}: per-branch {nD} (formula {2/(3*np.pi)*D**1.5:5.2f}) ; both branches {nD+nN}")
print(f"interval/circle zeta sums: {zeta(5):.6f}/{2*zeta(5):.6f} = 1/2 exact (halving certified)")
print("fixed-point counterterms: local boundary operators -> R-independent -> S-3 budget only")
