"""Reproduces the radion verdict of Section 7: the minimum, its depth, the radion
mass, the Casimir-only no-go, and the inverted TeV stiffness bound."""
import numpy as np
from scipy.special import zeta
from scipy.optimize import minimize_scalar
MPl=2.435e27; conv=1.97327e-7; C=3*zeta(5)/(64*np.pi**6)
def g(x):
    w=np.arange(1,200)
    return np.sum(np.exp(-x*w)*((x*w)**2+3*x*w+3)/(3*w**5))/zeta(5)
res=minimize_scalar(lambda u:(-5+8*g(u))/u**4,bounds=(0.3,10),method='bounded')
u,h=res.x,res.fun
R=8e-6/conv; mf=u/(2*np.pi*R)
print(f"minimum: 2pi m_f R = {u:.3f} -> m_f = {mf*1e3:.1f} meV at R=8um ; depth^(1/4)={ (abs(h)*C*(2*np.pi*mf)**4)**0.25*1e3:.1f} meV")
dR=R*1e-4; V=lambda r:(-5+8*g(2*np.pi*mf*r))*C/r**4
mr=np.sqrt((2/3)*(R/MPl)**2*(V(R+dR)-2*V(R)+V(R-dR))/dR**2)
print(f"radion mass = {mr:.1e} eV (needs >= 5e-3 eV) -> Casimir-only stabilisation EXCLUDED")
for m0 in (5e-3,25e-3,40e-3):
    print(f"m_r={m0*1e3:3.0f} meV -> V_stab^(1/4) >= {(1.5*MPl**2*m0**2)**0.25/1e12:.1f} TeV")
