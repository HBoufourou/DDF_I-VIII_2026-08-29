"""Reproduces the factorised S-3 system of §7: V_eff = -C_net/2R^4 + V_GW + V_bare."""
import numpy as np
from scipy.special import zeta
from scipy.optimize import brentq
MPl=2.435e27; conv=1.97327e-7; Cnet=3*3*zeta(5)/(64*np.pi**6)/2
M=25e-3; R=8.2e-6/conv; M5=(MPl**2/(np.pi*R))**(1/3); v2=M5**3; Lam=(2.25e-3)**4
VC=lambda r_: -Cnet/r_**4
def VGW(r_,q):
    L=np.pi*r_; return M*((v2+v2/q**2)*np.cosh(M*L)-2*v2/q)/np.sinh(M*L)
d=lambda q:(VC(R*1.000001)+VGW(R*1.000001,q)-VC(R*0.999999)-VGW(R*0.999999,q))
q=brentq(d,1.5,200.)
Vb=Lam-(VC(R)+VGW(R,q))
print(f"solution: v0/v1={q:.1f} ; V_bare=-({abs(Vb)**0.25/1e12:.2f} TeV)^4 ; tuning {Lam/abs(Vb):.1e}")
