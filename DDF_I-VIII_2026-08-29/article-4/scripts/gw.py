"""Reproduces the constructed flat-space Goldberger-Wise stabilisation of §7."""
import numpy as np
MPl=2.435e27; conv=1.97327e-7; M=25e-3
R=8e-6/conv; Ls=np.pi*R; M5=(MPl**2/Ls)**(1/3)
print(f"boundary hierarchy for R*=8um: v0/v1 = e^(M piR) = {np.exp(M*Ls):.1f}")
v0sq=M5**3; v1sq=v0sq/np.exp(2*M*Ls)
V=lambda L:M*((v0sq+v1sq)*np.cosh(M*L)-2*np.sqrt(v0sq*v1sq))/np.sinh(M*L)
dL=Ls*1e-5; Vpp=(V(Ls+dL)-2*V(Ls)+V(Ls-dL))/dL**2
print(f"V_stab^(1/4) = {abs(V(Ls))**0.25/1e12:.1f} TeV ; m_r = {np.sqrt((2/3)*(R/MPl)**2*Vpp*np.pi**2)*1e3:.1f} meV (floor: 5)")
print(f"V_min = M(v0^2-v1^2) > 0 : tunable uplift contribution available")
