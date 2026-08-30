"""POINT 1 (le Graal) — LVS a 2 modules de Kahler : la hierarchie passe-t-elle,
et que vaut c_b via le dictionnaire ? K = -2 ln(V + xi/2), V = tb^1.5 - l ts^1.5,
W = W0 + A e^{-a Ts}. Axions au minimum ; S/structures complexes gelees lourdes (flux)."""
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import eigh
xi,l,W0,A,a=0.5,1.0,1.0,1.0,2*np.pi/3
def Kf(t):
    tb,ts=t; V=tb**1.5-l*ts**1.5
    return -2*np.log(V+xi/2) if V+xi/2>0 else 1e9
def Wf(ts): return W0-A*np.exp(-a*ts)      # signe standard LVS (axion aligne)
def dWs(ts): return a*A*np.exp(-a*ts)*0.5  # d/dT_s = (1/2) d/dtau_s
def VF(t):
    tb,ts=t
    if tb<2 or ts<0.5 or tb**1.5<=l*ts**1.5+0.1: return 1e3
    h=1e-5
    K0=Kf(t)
    g=np.array([(Kf([tb+h,ts])-Kf([tb-h,ts]))/(2*h),(Kf([tb,ts+h])-Kf([tb,ts-h]))/(2*h)])
    Kt=np.zeros((2,2))
    for i in range(2):
        for j in range(2):
            e1=np.zeros(2);e1[i]=h;e2=np.zeros(2);e2[j]=h
            Kt[i,j]=(Kf(t+e1+e2)-Kf(t+e1-e2)-Kf(t-e1+e2)+Kf(t-e1-e2))/(4*h*h)
    Kij=Kt/4                                  # K_{T_i Tbar_j} = (1/4) d2K/dtau_i dtau_j
    DW=np.array([0.5*g[0]*Wf(ts), dWs(ts)+0.5*g[1]*Wf(ts)])
    try: Kinv=np.linalg.inv(Kij)
    except: return 1e3
    return np.exp(K0)*(DW@Kinv@DW-3*Wf(ts)**2)
best=None
for tb0 in (50,200,1000,5000):
    for ts0 in (2,4,7):
        m=minimize(VF,[tb0,ts0],method='Nelder-Mead',options={'xatol':1e-9,'fatol':1e-22,'maxiter':4000})
        if m.fun<1e2 and (best is None or m.fun<best.fun): best=m
tb,ts=best.x; Vol=tb**1.5-l*ts**1.5
print(f"minimum LVS : tau_b = {tb:.1f}, tau_s = {ts:.3f} ; Volume = {Vol:.1f} ; V_min = {best.fun:.2e}")
# masses : Hessienne de VF + metrique cinetique (2 K_ij pour tau reels) -> eig generalise
h=1e-4*np.array([tb,ts])
H=np.zeros((2,2))
for i in range(2):
    for j in range(2):
        e1=np.zeros(2);e1[i]=h[i];e2=np.zeros(2);e2[j]=h[j]
        H[i,j]=(VF(best.x+e1+e2)-VF(best.x+e1-e2)-VF(best.x-e1+e2)+VF(best.x-e1-e2))/(4*h[i]*h[j])
Kt=np.zeros((2,2)); hh=1e-5
for i in range(2):
    for j in range(2):
        e1=np.zeros(2);e1[i]=hh;e2=np.zeros(2);e2[j]=hh
        Kt[i,j]=(Kf(best.x+e1+e2)-Kf(best.x+e1-e2)-Kf(best.x-e1+e2)+Kf(best.x-e1-e2))/(4*hh*hh)
ev,evec=eigh(H,Kt/2)                        # L_cin = (K_tt/2)(d tau)^2 conventions coherentes
ev=np.real(ev)
print(f"m2 (base canonique) : leger {ev[0]:.2e}, lourd {ev[1]:.2e} ; HIERARCHIE = {abs(ev[0]/ev[1]):.2e}  {'PASSE (<0,1)' if abs(ev[0]/ev[1])<0.1 else 'ECHOUE'}")
# etat leger en fraction tau_b (canonique)
vl=evec[:,0]; nrm=np.sqrt(vl@ (Kt/2) @vl); vl=vl/nrm
fb=(vl[0]**2*(Kt[0,0]/2))/(vl@(Kt/2)@vl)
print(f"etat leger : fraction tau_b (grand volume) = {fb:.1%}  -> quasi pur volume global")
print()
print("=== LE DICTIONNAIRE APPLIQUE (l'etat leger ~ direction ln V6) ===")
a1=1/(4*np.sqrt(7)); b1=-np.sqrt(7)/4; b2=-1/np.sqrt(35)
gT=np.array([0.,5*b2,5*a1+b1]); nV=gT/np.linalg.norm(gT)
w=np.array([5*np.sqrt(2)/4,np.sqrt(10/7),9/(2*np.sqrt(14))])
cb=abs(nV@w); th=np.degrees(np.arccos(min(1,cb/np.linalg.norm(w))))
print(f"  c_b(LVS vanille) = |n_V6 . w| = {cb:.4f} = 7 sqrt6/12 (exact, la projection volume-global)")
print(f"  angle theta = {th:.1f} deg ; contenu radion = {nV[2]**2:.1%}")
print(f"  fenetre corpus [0,30 ; 0,61] : {'DEDANS' if 0.30<=cb<=0.61 else 'HORS FENETRE'} ; fenetre E-11 [73,4 ; 81,9] : {'DEDANS' if 73.4<=th<=81.9 else 'HORS'}")
