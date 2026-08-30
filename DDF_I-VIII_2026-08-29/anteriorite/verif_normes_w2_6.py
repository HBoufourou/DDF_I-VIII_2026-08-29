"""Anteriorite |w|^2 = 6 -- script de verification (23/08/2026).
Reproduit : (A) les normes du corpus ; (B) le transport de l'invariant Delta
de Lu-Pope ; (C) la dependance en p ; (D) la decomposition sur le plan
(tau, rho) de Wrase-Zagermann. Zero entree ajustable."""
import numpy as np
a1=1/(4*np.sqrt(7)); b1=-np.sqrt(7)/4; a2=np.sqrt(5/28); b2=-1/np.sqrt(35)
S2=np.sqrt(2)
w=np.array([5*S2/4, S2*(4*a2+5*b2), S2*9*a1])
v=np.array([5*S2/2, S2*(4*a2+5*b2), S2*(9*a1+b1)])
assert abs(w@w-6)<1e-12 and abs(v@v-14)<1e-12 and abs(w@v-8)<1e-12
D10=(5/2)**2+2*(-1)*9/8                      # Delta(F0) en 10D
c2=D10-2*(-1)*3/2                            # transporte en 4D (d=-1, dt=3)
assert abs(D10-4)<1e-12 and abs(2*c2-(v@v))<1e-12
for p in range(3,9):
    x=p-3
    wp=np.array([S2*x/4, S2*(4*a2+x*b2), S2*(p+1)*a1])
    assert abs(wp@wp-(x**2/5-x+6))<1e-12
g6=np.array([0,5*b2,5*a1+b1]); gt=np.array([-S2,5*b2/2,(5*a1+b1)/2]); gr=g6/3
e1=gt/np.linalg.norm(gt); t=gr-(gr@e1)*e1; e2=t/np.linalg.norm(t)
pr=(w@e1)*e1+(w@e2)*e2
assert abs(pr@pr-31/6)<1e-12 and abs(w@w-pr@pr-5/6)<1e-12
print("OK : |w|^2=6 ; |v|^2=14=2x(Delta=4 transporte) ; |w(p)|^2=(p-3)^2/5-(p-3)+6")
print("     6 atteint pour p=3 et p=8 seulement ; projection WZ = 31/6 ; residu = 5/6")
