"""POINT 2 — LE SCAN : le theta = 79,6 deg du premier vide est-il generique ou chanceux ?
Balayage des parametres du racetrack + criteres du CDD (hierarchie ET contenu radion)."""
import numpy as np
from scipy.optimize import minimize
a1=1/(4*np.sqrt(7)); b1=-np.sqrt(7)/4; b2=-1/np.sqrt(35)
gS=np.array([-np.sqrt(2),5*b2,5*a1]); gT=np.array([0.,5*b2,5*a1+b1])
nS=gS/np.linalg.norm(gS); t=gT-(gT@nS)*nS; nT=t/np.linalg.norm(t)
w=np.array([5*np.sqrt(2)/4,np.sqrt(10/7),9/(2*np.sqrt(14))]); wh=w/np.linalg.norm(w)
rng=np.random.default_rng(20260805)
def solve(A_,a_,B_,b_,W0):
    W=lambda s,t: W0+A_*np.exp(-a_*t)+B_*np.exp(-b_*s)
    def V(x):
        s,t=x
        if s<=0.05 or t<=0.05: return 1e3
        K=-np.log(2*s)-3*np.log(2*t)
        DS=-b_*B_*np.exp(-b_*s)-W(s,t)/(2*s); DT=-a_*A_*np.exp(-a_*t)-3*W(s,t)/(2*t)
        return np.exp(K)*((2*s)**2*DS**2+(2*t)**2/3*DT**2-3*W(s,t)**2)
    best=None
    for s0,t0 in ((5,5),(10,10),(20,15),(30,25)):
        m=minimize(V,[s0,t0],method='Nelder-Mead',options={'xatol':1e-11,'fatol':1e-18,'maxiter':3000})
        if best is None or m.fun<best.fun: best=m
    s0,t0=best.x
    if s0<=0.2 or t0<=0.2: return None
    h=1e-5; Vc=lambda y: V([np.exp(np.sqrt(2)*y[0])/2, np.exp(y[1]*np.sqrt(2/3))/2])
    y0=np.array([np.log(2*s0)/np.sqrt(2),np.sqrt(3/2)*np.log(2*t0)])
    H=np.zeros((2,2))
    for i in range(2):
        for j in range(2):
            e1=np.zeros(2);e1[i]=h;e2=np.zeros(2);e2[j]=h
            H[i,j]=(Vc(y0+e1+e2)-Vc(y0+e1-e2)-Vc(y0-e1+e2)+Vc(y0-e1-e2))/(4*h*h)
    ev,evec=np.linalg.eigh(H)
    if ev[0]<=0: return None
    return ev[0]/ev[1], evec[:,0]
rows=[]
for _ in range(400):
    A_=rng.uniform(0.5,2); B_=rng.uniform(0.5,2)
    N1,N2=rng.integers(8,16),rng.integers(8,16)
    if N1==N2: continue
    a_=2*np.pi/N1; b_=2*np.pi/N2; W0=-10**rng.uniform(-4,-2)
    r=solve(A_,a_,B_,b_,W0)
    if r is None: continue
    hier,e=r
    n=e[0]*nS+e[1]*nT
    th=np.degrees(np.arccos(min(1,abs(n@wh)))); cb=abs(n@w); rad=n[2]**2
    rows.append((hier,th,cb,rad))
R=np.array(rows)
print(f"vides stabilises : {len(R)}")
print(f"  theta : mediane {np.median(R[:,1]):.1f} deg ; quartiles [{np.percentile(R[:,1],25):.1f}, {np.percentile(R[:,1],75):.1f}]")
print(f"  c_b   : mediane {np.median(R[:,2]):.3f} ; quartiles [{np.percentile(R[:,2],25):.3f}, {np.percentile(R[:,2],75):.3f}]")
print(f"  contenu radion : mediane {np.median(R[:,3]):.1%}")
print(f"  hierarchie : mediane {np.median(R[:,0]):.3f} ; min {R[:,0].min():.3f}")
inw=(R[:,1]>=73.4)&(R[:,1]<=81.9)
print(f"\n  fraction dans la fenetre E-11 : {inw.mean():.1%}")
sel=inw&(R[:,0]<0.1)
print(f"  fraction qui passe AUSSI la hierarchie (<0,1) : {sel.mean():.1%}  ({sel.sum()} vides)")
print(f"\n  => le theta ~ 79,6 du premier point est-il generique ? "
      f"{'OUI (la famille est concentree la)' if abs(np.median(R[:,1])-79.6)<6 else 'NON : c est un point parmi une large dispersion'}")
