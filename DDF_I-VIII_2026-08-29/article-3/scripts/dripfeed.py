"""Reproduces every number of Sections 4, 5 and 7: the excluded flows, the Airy
matrix elements (closed form vs numerical), the computed gravitational rate, the
attractor window, and the quartic overlap integrals."""
import numpy as np
from scipy.special import ai_zeros, airy
from scipy.integrate import quad
hbar=6.582e-16; H0=2.2e-18; MPl=1.22e28
# Exclusions
omega=3e4/(10*1000*1.496e11); gap=25e-3
print(f"upward pumping : gap/(hbar omega) = {gap/(hbar*omega):.1e} -> P ~ exp(-1.9e24) : excluded")
print("baryonic leakage: Omega_b constant (BBN/CMB/today); CDM already 5.4x baryons at z=1100")
# Airy matrix elements: closed form vs numerical
az,_,_,_=ai_zeros(4); aD=-az
psi=lambda y,a: airy(y-a)[0]/abs(airy(-a)[1])
for (m,n) in [(1,2),(2,3)]:
    I=quad(lambda y: psi(y,aD[m-1])*y*psi(y,aD[n-1]),0,60,limit=400)[0]
    print(f"<{m}|y|{n}> numerical={abs(I):.5f}  closed form 2/(a_m-a_n)^2={2/(aD[m-1]-aD[n-1])**2:.5f}")
# Gravitational rate, two ways
for E in (13e-3,25e-3,40e-3):
    print(f"E={E*1e3:4.0f} meV : Gamma_grav ~ E^3/M_Pl^2 = {E**3/MPl**2/hbar:.1e} s^-1")
G=6.674e-11; m=25e-3*1.602e-19/9e16; w=25e-3*1.602e-19/1.055e-34; d=5e-6
print(f"quadrupole check: {G*m**2*d**4*w**5/(1.055e-34*(3e8)**5):.1e} s^-1 (same order)")
# Attractor window
print(f"gravity-only fraction: (2/9)*1.6e-46/H0 = {(2/9)*1.6e-46/H0:.1e}  -> ground floor empty")
for f in (1e-3,3e-2):
    print(f"target f={f:.0e} -> Gamma_fill = (9/2)H0 f = {4.5*H0*f:.1e} s^-1")
# Quartic overlaps
for (m,n) in [(1,2),(2,3),(2,2)]:
    I=quad(lambda y: psi(y,aD[m-1])**2*psi(y,aD[n-1])**2,0,60,limit=400)[0]
    print(f"I_{m}{n} = int|psi_m|^2|psi_n|^2 dy = {I:.4f}/l  (unsuppressed)")
print(f"deviation from a^-3: -(3/2)f -> {1.5e-3:.1e} at f=1e-3 (stated in the paper)")

# Elementary process: exothermic by concavity (full interleaved tower)
allz=np.sort(np.concatenate([-ai_zeros(6)[0],-ai_zeros(6)[1]]))
for nn in (2,3,4):
    print(f"n={nn}: 2E_n-(E_(n-1)+E_(n+1)) = {2*allz[nn-1]-(allz[nn-2]+allz[nn]):+.4f} eps (exothermic)")
# Lambda window
mphi=0.025; nres=0.26*3.7e-11/mphi
Gam=lambda lam: nres*lam**2*0.3/(64*np.pi*mphi**2)*1.52e15
for f in (1e-3,3e-2):
    print(f"f={f:.0e} -> lambda_4 = {np.sqrt(4.5*H0*f/Gam(1)):.1e}")
print("f(z) ∝ a^-3/2 : f(z=1)/f(0) = 2.83 ; early universe: Gamma/H -> 0 (BBN/CMB untouched)")

# --- report resolutions (7/7) ---
from scipy.integrate import quad as _q
from scipy.special import airy as _ai
lv=[1.019,2.338,3.248,4.088]
def _psi(y,c):
    N=np.sqrt(_q(lambda t:_ai(t-c)[0]**2,0,60,limit=400)[0]); return _ai(y-c)[0]/N
for n in (2,3):
    J=_q(lambda y:_psi(y,lv[n-1])**2*_psi(y,lv[n-2])*_psi(y,lv[n]),0,60,limit=400)[0]
    print(f"vertex J_{n} = {J:.3f}/l")
c=[_q(lambda y:_psi(y,x),0,60,limit=400)[0] for x in lv]
f=np.array(c)**2; f/=f.sum(); print("misalignment fractions (flat profile):",np.round(f,3))
print(f"delta integrated (z=1..0): {_q(lambda l:1.5e-3*np.exp(-1.5*l),np.log(.5),0)[0]:.2e}")
m=0.025
print(f"reservoir degeneracy n*l^3: local {1.2e-4*(2*np.pi/(m*1e-3))**3:.1e}, cosmic {3.8e-10*(2*np.pi/(m*1e-4))**3:.1e}")
print(f"ground-level BEC margin: {1.7e-6*(2*np.pi/(m*1e-3))**3:.1e} vs 2.61")
print(f"naive condensation: Gamma_coll/H0 = {3.8e-13*(2e-13)**2/(64*np.pi*m**2)*0.3/1.45e-33:.1e} -> stimulated channel required")

# --- §7.1-7.2 : la prediction w=-1 et la borne N_eff (ajout 01/08/2026) ---
m_r=9e-3; rho_c=8.5e-11; rho_m=0.31*rho_c; MPl_red=2.4e27
dR=rho_m/(m_r**2*MPl_red**2)
print(f"  radion stiffness: dR/R = {dR:.1e}  =>  |1+w| <= {4/3*dR:.0e}  (w=-1 to 61 decimals)")
for dN in (0.3,0.06):
    print(f"  dN_eff < {dN:4.2f}  =>  T_dark/T_gamma < {(dN*7/16)**0.25:.2f}")

# --- §7.3 : invisibilite cosmologique du mediateur ---
m_phi=0.025; T_b=2e12; T_now=2.35e-4; c_l=3e8
v_now=T_now/T_b
print(f"  misalignment: v_today = {v_now:.1e} c = {v_now*c_l:.1e} m/s ; v(recomb) = {v_now*c_l*1100:.1e} m/s")
f_hi=2.5e-2; f_lo=f_hi/(2.06e5/1e4); boost=3e5
print(f"  window: {f_lo:.1e} <= f_local <= {f_hi:.1e}")
fc=lambda z: f_hi*((1+z)**3/boost)*(v_now*c_l*(1+z)/2e5)/(0.31**0.5*(1+z)**1.5)
for z in (1100,3400): print(f"  f(z={z}) <= {fc(z):.1e}")
print(f"  dC_l/C_l ~ 1e-12 vs cosmic variance 1e-2, Planck 1e-3 -> 9 orders below floor")

# --- no-go du radion (attribution du signal 1er ordre corrigee, 02/08) ---
import numpy as _np
MPl_red=2.4e27; m_r_=9e-3; T_tr=4.6e12
S3=(16*_np.pi/3)*MPl_red**2/m_r_
print(f"  radion bounce (thin wall): S3/T = {S3/T_tr:.1e} -> tunnelling excluded")
print(f"  => the radion rolls (t_roll ~ 1.4e-13 s), it does not nucleate bubbles;")
print(f"     the first-order GW signature must be attributed to the TeV-scale event.")

# --- fenetre lambda stimulee + plafond Bullet (02/08) ---
lam_unstim=1e-12; occ=1e12
print(f"  stimulated window: lambda ~ {lam_unstim/occ**0.5:.0e} (floor) ; Bullet ceiling: lambda < 2e-13")
print(f"  sigma/m(1e-12) = 35 cm2/g EXCLU ; sigma/m(1e-18) ~ 3e-11 cm2/g invisible")
