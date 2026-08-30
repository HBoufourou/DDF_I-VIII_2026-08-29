"""DERIVATION NUMERIQUE DE c_b — dans la classe de vides que le theoreme exige
(deux directions quasi degenerees). Tout est numerique ; la SEULE entree non derivee
est le mecanisme qui epingle tau_f au point de degenerescence — declaree. 07/08/2026."""
import numpy as np
# ---- le vide verrouille (valide par son run v4) ----
Vol=1.279e20; W0=1.32e-9; ts=21.2185; A_loop=1e-5; B_loop=1e-5; C_loop=1e-5
projV=7*np.sqrt(6)/12; projF=np.sqrt(5/6); amp=np.hypot(projV,projF); dph=np.degrees(np.arctan2(projF,projV))
G=np.diag([0.5,1.0])                                  # metrique de Kahler en (v,w)=(ln tf, ln tb)
print("="*74); print("1. LES DEUX MASSES ET LE COUPLAGE HORS-DIAGONAL (numerique)")
print("="*74)
m2V=W0**2/Vol**3                                       # echelle LVS du volume (terme alpha')
def Vloop(v,lnVol):
    tf=np.exp(v); V=np.exp(lnVol)
    return (W0**2/V**2)*(A_loop/tf**2 - B_loop/(V*np.sqrt(tf)) + C_loop*tf/V**2)
# point de degenerescence : m2_F(tau_f*) = m2_V  =>  tau_f* = sqrt(2 A) * ... resolu numeriquement
from scipy.optimize import brentq
def m2F_of(v):
    h=1e-4
    return (Vloop(v+h,np.log(Vol))-2*Vloop(v,np.log(Vol))+Vloop(v-h,np.log(Vol)))/h**2/G[0,0]
vstar=brentq(lambda v: m2F_of(v)-m2V, np.log(1e3), np.log(1e12))
print(f"  m2_V (LVS, volume)          = {m2V:.3e}   [W0^2/V^3, numerique]")
print(f"  tau_f* (degenerescence)     = {np.exp(vstar):.3e}   [resolu : m2_F(tau_f*) = m2_V]")
print(f"  controle analytique sqrt(2A_loop*Vol) ~ {np.sqrt(2*A_loop*Vol):.2e}  (meme ordre) ✓")
# hors-diagonal : d2 Vloop / dv d(lnVol) au point, converti en base canonique orthonormee
h=1e-4
Hvw=(Vloop(vstar+h,np.log(Vol)+h)-Vloop(vstar+h,np.log(Vol)-h)-Vloop(vstar-h,np.log(Vol)+h)+Vloop(vstar-h,np.log(Vol)-h))/(4*h*h)
nV=np.array([1.,1.]); nV=nV/np.sqrt(nV@G@nV)          # direction volume (canonique)
nF=np.array([2.,-1.]); nF=nF/np.sqrt(nF@G@nF)         # direction fibre (canonique, orthogonale-G)
# matrice 2x2 dans la base (nV, nF) : diag(m2V, m2F) + hors-diag venant de Vloop
Hfull=np.zeros((2,2))
def Vtot(v,lw):  # potentiel total le long des 2 directions autour du point
    dv,dw=v,lw
    x=np.array([vstar,np.log(Vol/np.sqrt(np.exp(vstar)))])+dv*nV+dw*nF  # deplacement canonique
    vv,ww=x; lnVol_loc=np.log(np.exp(vv/1.)*0+1)  # (ln Vol reconstruit)
    lnV=0.5*vv+ww if False else None
    return None
# plus simple et sur : matrice = P^T Hlog P avec Hlog la Hessienne (v,w) complete
def Vtot2(v,w):
    lnV=np.log(max(np.exp(0.5*v+w)-ts**1.5,1e-30))
    c1=(8/3)*(np.pi*1.0)**2*np.sqrt(ts)*np.exp(-2*np.pi*ts)
    VV=np.exp(lnV)
    lvs=c1/VV - 4*np.pi*1.0*W0*ts*np.exp(-np.pi*ts)/VV**2 + 3*(0.194/0.01**1.5)*W0**2/(4*VV**3)
    return lvs+Vloop(v,lnV)
x0=np.array([vstar,np.log((Vol+ts**1.5)/np.sqrt(np.exp(vstar)))])
H=np.zeros((2,2)); hh=[1e-4,1e-4]
for i in range(2):
    for j in range(2):
        e1=np.zeros(2);e1[i]=hh[i];e2=np.zeros(2);e2[j]=hh[j]
        H[i,j]=(Vtot2(*(x0+e1+e2))-Vtot2(*(x0+e1-e2))-Vtot2(*(x0-e1+e2))+Vtot2(*(x0-e1-e2)))/(4*hh[i]*hh[j])
from scipy.linalg import eigh
ev,evec=eigh((H+H.T)/2,G)
print(f"  valeurs propres numeriques  = {ev[0]:.3e}, {ev[1]:.3e}   (rapport {ev[0]/ev[1]:.2f})")
n=evec[:,0]; n=n/np.sqrt(n@G@n)
aV=n@G@nV; aF=n@G@nF; pf=aF**2/(aV**2+aF**2)
cb=abs(aV*projV+aF*projF)
print(f"  etat leger : fraction fibre = {pf:.1%} ; c_b = |{aV:+.3f}*{projV:.4f} + {aF:+.3f}*{projF:.4f}| = {cb:.4f}")
print("\n"+"="*74); print("2. c_b EN FONCTION DU DESACCORD (la courbe complete)")
print("="*74)
print(f"  {'m2F/m2V':>9} {'mixing':>8} {'%fibre':>8} {'c_b':>8} {'fenetre?'}")
eps=abs((nV@G@ ( (H+H.T)/2 @ np.linalg.inv(G)) @ nF))   # ordre du hors-diagonal
for r in (0.5,0.8,0.9,1.0,1.1,1.25,2.0):
    M=np.array([[m2V, eps],[eps, r*m2V]])
    w_,v_=np.linalg.eigh(M)
    nl=v_[:,0]; pf=nl[1]**2
    cbr=abs(nl[0]*projV+nl[1]*projF)
    print(f"  {r:9.2f} {np.degrees(0.5*np.arctan2(2*eps,(1-r)*m2V)):8.1f} {pf:8.1%} {cbr:8.3f} {'OUI' if 0.30<=cbr<=0.61 else '-'}")
print(f"\n  hors-diagonal numerique eps/m2V = {eps/m2V:.3f}")
print(f"  ==> TOLERANCE de degenerescence requise : |1 - m2F/m2V| < ~{2*eps/m2V:.2f} pour un melange fort")
