"""EXP-B v1.1 — J1'/J2' : capaciteur pair, deviation du profil ET cadran des RAPPORTS.
Remplace le bloc numerique casse (100 %) de J1prime.py v2. Boufourou, 03/08/2026."""
import numpy as np
from scipy.linalg import eigh_tridiagonal

hbar_c=197.3      # meV.um
ell=5.802         # um (echelle d'Airy du corpus)
u=25.76/ell       # intervalle en unites d'Airy = 4.44 (le u du corpus)

# ---------- 1. deviation du PROFIL (ce que son script devait faire) ----------
def profil_dev(m_meV):
    mu=m_meV*ell/hbar_c              # m_sigma * ell, sans dimension
    x=np.linspace(0,u,2000)
    V_lin=x
    V_m=np.sinh(mu*x)/mu if mu>0 else x   # pente unitaire en 0 (meme k)
    return np.max(np.abs(V_m-V_lin))/np.max(V_lin)

# ---------- 2. le VRAI cadran : les rapports de la tour ENTRELACEE ----------
def tour(m_meV,N=3000):
    mu=m_meV*ell/hbar_c
    h=u/N; x=np.linspace(h,u-h,N-1)
    V=np.sinh(mu*x)/mu if mu>0 else x
    d=1/h**2+V; e=-1/(2*h**2)*np.ones(N-2)
    # branche impaire : Dirichlet aux deux bords
    Eo=eigh_tridiagonal(d,e,select='i',select_range=(0,3))[0]
    # branche paire : Neumann en 0 (point miroir), Dirichlet en u
    d2=d.copy(); d2[0]=1/(2*h**2)+V[0]   # sigma'(0)=0 via ghost
    Ee=eigh_tridiagonal(d2,e,select='i',select_range=(0,3))[0]
    E=np.sort(np.concatenate([Ee,Eo]))[:4]
    return E/E[0]

print("m_sigma  m.piR   dev profil (num vs (m.piR)^2/6)   rapports tour        ecart max rapports")
r0=tour(0)
print(f"  0.0    0.000        0.000 %  /  0.000 %         {np.round(r0,3)}   —   (corpus : 1, 2.29, 3.19, 4.01)")
for m in (0.5,1.0,1.5,2.0,3.0,5.0):
    mu_piR=m*25.76/hbar_c
    dev=profil_dev(m); anal=mu_piR**2/6
    r=tour(m); shift=np.max(np.abs(r-r0)/r0)
    print(f"  {m:3.1f}    {mu_piR:.3f}        {dev*100:6.3f} %  /  {anal*100:6.3f} %         {np.round(r,3)}       {shift*100:6.3f} %")
# seuil : rapports a 1 %
from scipy.optimize import brentq
f=lambda m: np.max(np.abs(tour(m)-r0)/r0)-0.01
m1=brentq(f,0.5,8)
print(f"\n==> CADRAN (sur les RAPPORTS, la quantite mesurable) : rapports a 1 %  =>  m_sigma < {m1:.1f} meV")
print(f"    (le seuil sur le PROFIL donnait 1.9 meV ; les rapports sont moins sensibles — le vrai seuil est {m1:.1f})")
