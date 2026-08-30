"""Article II — verification de chaque nombre du texte (C3 exact + route gravitationnelle + w1)."""
import numpy as np
from scipy.linalg import eigh_tridiagonal
hbar_c=197.3  # meV.um
piR=25.76; m_phi=24.0  # meV
F_mevum=24.0/5.802     # meV/um (pente : epsilon/ell, corpus)
ell=((hbar_c)**2/(2*m_phi*F_mevum))**(1/3)
u=piR/ell
print(f"ell = {ell:.3f} um (texte 5,81) ; u = piR/ell = {u:.3f} (texte 4,43)")
F_GeV2=8.161e-22
print(f"F = {F_GeV2:.2e} GeV^2 ; sqrt(F) = {np.sqrt(F_GeV2)*1e12:.1f} meV (texte 29)")
# route gravitationnelle : T^(1/4) requis
M5=3.568e8; mphi_GeV=24e-12  # convention piR (orbifold)
T=6*F_GeV2*M5**3/mphi_GeV
print(f"route grav : T^(1/4) = (6 F M5^3 / m_phi)^(1/4) = {T**0.25/1e3:.1f} TeV (texte 9,8)")
# tour entrelacee + w1
def tour(u,N=4000):
    h=u/N; x=np.linspace(h,u-h,N-1); V=x
    d=1/h**2+V; e=-1/(2*h**2)*np.ones(N-2)
    Eo=eigh_tridiagonal(d,e,select='i',select_range=(0,1))[0]
    d2=d.copy(); d2[0]=1/(2*h**2)+V[0]
    Ee,ve=eigh_tridiagonal(d2,e,select='i',select_range=(0,1))
    E=np.sort(np.concatenate([Ee,Eo]))[:4]
    psi=ve[:,0]; psif=np.concatenate([[psi[0]],psi]); xf=np.concatenate([[0],x])
    psif/=np.sqrt(np.trapezoid(psif**2,xf))
    return E/E[0], u*psif[0]**2
r,w1=tour(4.434)
print(f"rapports a u=4.434 : {np.round(r,3)} (texte 1:2.293:3.188:4.028)")
print(f"w1 = u|psi1(0)|^2 = {w1:.2f} (texte 5,5) ; alpha_b = {np.sqrt(2*w1):.2f} c_b (texte 3,3)")
for a in (1.0,2.0,10,512): print(f"  alpha_b requis {a:5} -> c_b = {a/np.sqrt(2*w1):.2f}")
# dial m_sigma : deja etabli 4,3 meV (J1J2_capaciteur.py) ; verif rapide 4e rapport vs u
for uu in (4.2,4.434,4.7):
    rr,_=tour(uu); print(f"  u={uu:.3f} : 4e rapport = {rr[3]:.3f} (sensibilite du dial 2)")
# tensions
J0max=4e-29; g_min=2*F_GeV2/J0max
print(f"J0max = {J0max:.0e} ; g_min = {g_min:.1e} GeV^(1/2) = {g_min/np.sqrt(M5):.0f} x sqrt(M5) (texte ~2e3)")
print(f"m_r(sigma) = (30 meV)^2/MPl = {(30e-12)**2/2.4e18:.0e} GeV = {(30e-12)**2/2.4e18*1e9:.0e} eV (texte 4e-31)")
