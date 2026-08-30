# ══════════════════════════════════════════════════════════════════════════════
#  ROUTE 1 — BINNING DE |p| PAR HIÉRARCHIE LVS SUR LES 20 POLYTOPES
#  CELLULE COLAB AUTONOME. Lit les κ dans les pickles CQSV : rien à envoyer.
#  Le cône de Kähler explicite n'est PAS nécessaire — l'admissibilité est
#  testée par la positivité de la métrique G, équivalente en pratique.
# ══════════════════════════════════════════════════════════════════════════════
import os, glob, gzip, pickle, itertools, subprocess, warnings, json
import numpy as np
warnings.filterwarnings("ignore")

# --- LES 20 POLYTOPES (scan CQSV définitif du 14/08) -------------------------
CIBLE_81 = [(3181,5),(3183,3),(3214,0),(3215,22)]                    # (h12,chi)=(81,-152)
NEUFS    = [(156,2),(373,0),(1093,5),(1099,3),(1128,22),(1234,1),(1243,0),
            (1839,2),(1950,0),(1998,0),(2439,0),(2899,0),(3662,0),
            (4255,5),(4855,5),(4856,5)]                              # Hodge différents
TOUS = CIBLE_81 + NEUFS
W_GEOM, FEN, VT, NTIR = 1.6956, (0.30,0.61), 1e6, 200000
REPO = "CY_Orientifold_database"

if not os.path.isdir(REPO):
    subprocess.run(["git","clone","--depth","1",
        "https://github.com/AndreasSchachner/CY_Orientifold_database.git"], check=False)
FOLDER = next((c for pat in ("**/h11_5","**/h11=5") for c in glob.glob(pat, recursive=True)
               if os.path.isdir(c) and glob.glob(os.path.join(c,"POLY_*.p"))), None)
print(f"dossier : {FOLDER}\n")
L = lambda x: x.tolist() if hasattr(x,"tolist") else list(x)

def structure(P,T):
    """K3, dP diagonaux, Wilson — lus dans DIVS/DPEZ/SCHERN (déjà calculés par les auteurs)"""
    types={}
    for row,c2 in zip(L(P['DIVS'][T]), L(P['SCHERN'][T])):
        r=L(row); i=int(r[0]); h=[int(v) for v in L(r[1])]; chi_h=float(r[4])
        if h[:4]==[1,0,1,20] and int(c2)==24: types[i]="K3"
        elif chi_h==0.0:                      types[i]="Wilson"
        elif h[1]==0 and h[2]==0:             types[i]=f"dP{h[3]-1}"
        else:                                  types[i]="autre"
    ddP=[int(L(e)[0]) for e in L(P['DPEZ'][T]) if len(L(e))>2 and int(L(e)[2])==1]
    K3=[i for i,t in types.items() if t=="K3"]
    return (K3, ddP, types) if (K3 and len(ddP)>=2) else (None,None,None)

def kappa(P,T,basis):
    k=np.zeros((5,5,5))
    for key,val in dict(P['INTNUMS'][T]).items():
        t=tuple(int(x) for x in (key if isinstance(key,tuple)
                else str(key).strip("()").replace(" ","").split(",")))
        if all(x in basis for x in t):
            a=tuple(basis.index(x) for x in t)
            for p in set(itertools.permutations(a)): k[p]=float(val)
    return k

def mesure(kap,t,iK,f1,f2):
    """renvoie (|p|, tau_dP/V^(2/3)) si le point est admissible, sinon None"""
    V=np.einsum('ijk,i,j,k',kap,t,t,t)/6.
    if V<=0: return None
    M=np.einsum('ijk,k',kap,t)
    if abs(np.linalg.det(M))<1e-25: return None
    tau=np.einsum('ijk,j,k',kap,t,t)/2.
    if np.any(tau<=0): return None
    G=-np.linalg.inv(M)/V+np.outer(t,t)/(2*V**2); G=0.5*(G+G.T)
    ev=np.linalg.eigvalsh(G)
    if ev.min()<=0 or ev.max()/ev.min()>1e14: return None    # <- admissibilité
    S=np.linalg.cholesky(G)
    A=np.array([S.T@f for f in [t/(2*V),np.eye(5)[f1],np.eye(5)[f2]]])
    Q,_=np.linalg.qr(A.T); u,s,_=np.linalg.svd(np.eye(5)-Q@Q.T)
    B=[np.linalg.solve(S.T,u[:,j]) for j in range(5) if s[j]>1e-8][:2]
    if len(B)<2: return None
    B=[x/np.sqrt(x@G@x) for x in B]
    cov=np.eye(5)[iK]/tau[iK]; nr=np.sqrt(cov@np.linalg.inv(G)@cov)
    p=float(np.hypot(*[W_GEOM*(cov@x)/nr for x in B]))
    if not np.isfinite(p) or p>W_GEOM*1.0001: return None    # <- PLAFOND
    return p, min(tau[f1],tau[f2])/V**(2/3.)

BINS=[(1e-1,1e9,"~1 (pas LVS)"),(1e-2,1e-1,"1e-2"),(1e-3,1e-2,"1e-3"),
      (1e-4,1e-3,"1e-4"),(1e-6,1e-4,"1e-6..1e-4"),(0,1e-6,"<1e-6 (LVS)")]
res={}
for pid,T in TOUS:
    fp=glob.glob(os.path.join(FOLDER,f"POLY_{pid}_*.p"))
    if not fp: print(f"POLY {pid} : absent"); continue
    with gzip.open(fp[0],"rb") as f: P=pickle.load(f,encoding="latin1")
    K3,ddP,types=structure(P,T)
    if K3 is None: print(f"POLY {pid} tri {T} : structure absente"); continue
    D=None
    for e1,e2 in itertools.combinations([x for x in range(1,10) if x not in [K3[0]]+ddP[:2]],2):
        basis=sorted([K3[0],ddP[0],ddP[1],e1,e2])
        kap=kappa(P,T,basis)
        if np.count_nonzero(kap)<5: continue
        iK,f1,f2=basis.index(K3[0]),basis.index(ddP[0]),basis.index(ddP[1])
        rng=np.random.default_rng(0); dat=[]
        for _ in range(NTIR):
            t=np.exp(rng.uniform(-5,4,5))
            Vc=np.einsum('ijk,i,j,k',kap,t,t,t)/6.
            if Vc<=0: continue
            t*=(VT/Vc)**(1/3.)
            r=mesure(kap,t,iK,f1,f2)
            if r: dat.append(r)
        if len(dat)>=500 and (D is None or len(dat)>len(D[1])): D=(basis,np.array(dat))
    if D is None: print(f"POLY {pid} tri {T} : aucune base exploitable"); continue
    basis,dat=D
    print(f"\n=== POLY {pid} tri {T}  h12={int(P['h12'])}  base {basis}  "
          f"{len(dat)} points admissibles")
    print(f"  {'tau_dP/V^(2/3)':>16} {'n':>7} {'|p| median':>11} {'quartiles':>18} {'fenetre':>8}")
    row={}
    for lo,hi,lab in BINS:
        m=(dat[:,1]>=lo)&(dat[:,1]<hi)
        if m.sum()<20: continue
        v=dat[m,0]; q1,q3=np.percentile(v,[25,75])
        row[lab]=dict(n=int(m.sum()),med=float(np.median(v)),q1=float(q1),q3=float(q3),
                      frac=float(np.mean((v>=FEN[0])&(v<=FEN[1]))))
        print(f"  {lab:>16} {m.sum():7d} {np.median(v):11.4f} "
              f"{f'[{q1:.3f} ; {q3:.3f}]':>18} {row[lab]['frac']:7.1%}")
    res[f"{pid}_{T}"]=dict(h12=int(P['h12']),chi=int(P['chi']),basis=basis,bins=row)

json.dump(res,open("route1_binning.json","w"),indent=1)
print("\n" + "="*70)
print("SYNTHESE — mediane de |p| dans le regime LVS (tau_dP/V^(2/3) < 1e-3)")
print("="*70)
print(f"  {'POLY':>7} {'h12':>5} {'median LVS':>11} {'dans la fenetre ?':>18}")
for k,v in res.items():
    b=v['bins'].get('<1e-6 (LVS)') or v['bins'].get('1e-6..1e-4') or v['bins'].get('1e-4')
    if not b: continue
    ok = "OUI" if FEN[0]<=b['med']<=FEN[1] else "non"
    print(f"  {k.split('_')[0]:>7} {v['h12']:5d} {b['med']:11.4f} {ok:>18}")
print("\n-> route1_binning.json")
try:
    from google.colab import files; files.download("route1_binning.json")
except Exception: pass
