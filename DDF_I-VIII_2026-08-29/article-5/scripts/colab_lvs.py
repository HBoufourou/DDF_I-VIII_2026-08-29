# ══════════════════════════════════════════════════════════════════════════════
#  A4 — |p| AU POINT LVS  ·  CELLULE COLAB AUTONOME  ·  RIEN À ENVOYER
#  On ne tire plus au hasard : on RÉSOUT les 3 conditions physiques et on
#  balaie les 2 directions plates restantes. C'est ce qui marchait sur 4686.
# ══════════════════════════════════════════════════════════════════════════════
import os, glob, gzip, pickle, itertools, subprocess, warnings
import numpy as np
from scipy.optimize import fsolve
warnings.filterwarnings("ignore")

CIBLES = [(3181,5),(3181,7),(3183,3),(3215,22),(1234,1),(156,2),(1950,0),(373,0)]
W_GEOM, FEN = 1.6956, (0.30, 0.61)
VT, RS, TW  = 1e6, 1e-3, 1.51   # rs relache : |p| varie de 0,09 % de 1e-3 a 8,6e-12
REPO = "CY_Orientifold_database"

if not os.path.isdir(REPO):
    subprocess.run(["git","clone","--depth","1",
        "https://github.com/AndreasSchachner/CY_Orientifold_database.git"], check=False)
FOLDER = next((c for pat in ("**/h11_5","**/h11=5") for c in glob.glob(pat, recursive=True)
               if os.path.isdir(c) and glob.glob(os.path.join(c,"POLY_*.p"))), None)
print(f"dossier : {FOLDER}\n")
L = lambda x: x.tolist() if hasattr(x,"tolist") else list(x)

def dtype(h):
    a,b,c,d = [int(v) for v in L(h)[:4]]
    if b==0 and c==1 and d==20: return "K3"
    if b==1 and c==0:           return "Wilson"
    if b==0 and c==0:           return f"dP{d-1}"
    return "autre"

def build(P,T):
    types={}
    for row,c2 in zip(L(P['DIVS'][T]), L(P['SCHERN'][T])):
        r=L(row); types[int(r[0])] = "K3" if (dtype(r[1])=="K3" and int(c2)==24) else dtype(r[1])
    ddP=[int(L(e)[0]) for e in L(P['DPEZ'][T]) if len(L(e))>2 and int(L(e)[2])==1]
    K3=[i for i,t in types.items() if t=="K3"]; W=[i for i,t in types.items() if t=="Wilson"]
    if not(K3 and len(ddP)>=2 and W): return None
    for extra in range(1,12):
        base=[K3[0],ddP[0],ddP[1],W[0]]
        if extra in base: continue
        base=base+[extra]; kap=np.zeros((5,5,5))
        for key,val in dict(P['INTNUMS'][T]).items():
            k=tuple(int(x) for x in (key if isinstance(key,tuple)
                    else str(key).strip("()").replace(" ","").split(",")))
            if all(x in base for x in k):
                a=tuple(base.index(x) for x in k)
                for p in set(itertools.permutations(a)): kap[p]=float(val)
        if np.count_nonzero(kap)>3: return kap, base
    return None

V_  = lambda kap,t: np.einsum('ijk,i,j,k',kap,t,t,t)/6.
TAU = lambda kap,t: np.einsum('ijk,j,k',kap,t,t)/2.

def solve_pt(kap, a, b, seed):
    """impose V, tau[1] (le petit dP) et tau[3] (le Wilson) ; a,b = les 2 libres"""
    def eqs(x):
        t=np.array([x[0],x[1],a,x[2],b]); tau=TAU(kap,t)
        return [V_(kap,t)/VT-1., tau[1]/(RS*VT**(2/3.))-1., tau[3]/(TW*VT**(2/3.))-1.]
    x=fsolve(eqs,seed,xtol=1e-14,maxfev=8000)
    t=np.array([x[0],x[1],a,x[2],b])
    return (t,x) if np.max(np.abs(eqs(x)))<1e-8 else (None,None)

def p_of(kap,t):
    V=V_(kap,t)
    if V<=0: return None
    M=np.einsum('ijk,k',kap,t)
    if abs(np.linalg.det(M))<1e-30: return None
    G=-np.linalg.inv(M)/V+np.outer(t,t)/(2*V**2); G=0.5*(G+G.T)
    tau=TAU(kap,t); ev=np.linalg.eigvalsh(G)
    if ev.min()<=0 or np.any(tau<=0) or ev.max()/ev.min()>1e14: return None
    S=np.linalg.cholesky(G)
    A=np.array([S.T@f for f in [t/(2*V),np.eye(5)[1],np.eye(5)[3]]])
    Q,_=np.linalg.qr(A.T); u,s,_=np.linalg.svd(np.eye(5)-Q@Q.T)
    B=[np.linalg.solve(S.T,u[:,j]) for j in range(5) if s[j]>1e-8][:2]
    if len(B)<2: return None
    B=[x/np.sqrt(x@G@x) for x in B]
    cov=np.eye(5)[0]/tau[0]; nr=np.sqrt(cov@np.linalg.inv(G)@cov)
    p=float(np.hypot(*[W_GEOM*(cov@x)/nr for x in B]))
    return p if (np.isfinite(p) and p<=W_GEOM*1.0001) else None

print(f"{'POLY':>6} {'tri':>4} {'h12':>4} {'base':>18} {'sol':>5} {'p_min':>7} "
      f"{'p_max':>7} {'p_med':>7} {'fenetre':>8}")
for pid,T in CIBLES:
    fp=glob.glob(os.path.join(FOLDER,f"POLY_{pid}_*.p"))
    if not fp: print(f"  {pid:>6} : absent"); continue
    with gzip.open(fp[0],"rb") as f: P=pickle.load(f,encoding="latin1")
    r=build(P,T)
    if r is None: print(f"  {pid:6d} {T:4d} : structure absente"); continue
    kap,base=r; vals=[]; S0=VT**(1/3.)
    seeds=[np.array([1.8*S0,0.55*S0,0.55*S0]), np.array([S0,S0,S0]),
           np.array([3*S0,0.2*S0,0.2*S0]), np.array([0.5*S0,2*S0,0.5*S0])]
    for a in np.linspace(0.05*S0,3.0*S0,24):
        for b in np.linspace(0.05*S0,3.0*S0,24):
            for sd in seeds:
                t,x=solve_pt(kap,a,b,sd)
                if t is None: continue
                v=p_of(kap,t)
                if v is not None:
                    vals.append(v)
                    if x is not None: seeds[0]=x
                    break
    if len(vals)<5:
        print(f"  {pid:6d} {T:4d} {int(P['h12']):4d} {str(base):>18} {len(vals):5d}  "
              f"(trop peu de solutions)"); continue
    v=np.array(vals); fr=float(np.mean((v>=FEN[0])&(v<=FEN[1])))
    print(f"  {pid:6d} {T:4d} {int(P['h12']):4d} {str(base):>18} {len(v):5d} "
          f"{v.min():7.4f} {v.max():7.4f} {np.median(v):7.4f} {fr:7.1%}")
print("\nColle ce tableau dans le chat.")
