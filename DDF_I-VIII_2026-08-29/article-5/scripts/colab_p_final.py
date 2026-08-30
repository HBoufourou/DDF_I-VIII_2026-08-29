# ══════════════════════════════════════════════════════════════════════════════
#  A4 — |p| DES CANDIDATS  ·  CELLULE COLAB AUTONOME  ·  AUCUN FICHIER À ENVOYER
#  Elle reclone la base (gratuit), n'ouvre QUE les 8 polytopes utiles, et
#  calcule |p| avec une ASSERTION qui rend le bug précédent impossible.
#  Sortie : un petit tableau TEXTE, à copier-coller dans le chat.
# ══════════════════════════════════════════════════════════════════════════════
import os, glob, gzip, pickle, itertools, subprocess, warnings
import numpy as np
warnings.filterwarnings("ignore")

CIBLES = [(3181,5),(3181,7),(3183,3),(3215,22),(1234,1),(156,2),(1950,0),(373,0)]
W_GEOM = 1.6956            # |w_geom| — PLAFOND MATHÉMATIQUE de |p|
FEN    = (0.30, 0.61)
REPO   = "CY_Orientifold_database"

if not os.path.isdir(REPO):
    subprocess.run(["git","clone","--depth","1",
        "https://github.com/AndreasSchachner/CY_Orientifold_database.git"], check=False)
FOLDER = next((c for pat in (f"**/h11_5", f"**/h11=5")
               for c in glob.glob(pat, recursive=True)
               if os.path.isdir(c) and glob.glob(os.path.join(c,"POLY_*.p"))), None)
print(f"dossier : {FOLDER}\n")

def L(x):
    try: return x.tolist()
    except Exception: return list(x)

def dtype(h):
    h00,h10,h20,h11 = [int(v) for v in L(h)[:4]]
    if h10==0 and h20==1 and h11==20: return "K3"
    if h10==1 and h20==0:             return "Wilson"
    if h10==0 and h20==0:             return f"dP{h11-1}"
    return "autre"

def build(P, T):
    """(kappa 5x5x5 FRAIS, index K3, index dP diagonal, index Wilson) — tout relu ici"""
    types = {}
    for row,c2 in zip(L(P['DIVS'][T]), L(P['SCHERN'][T])):
        r = L(row); i = int(r[0])
        types[i] = "K3" if (dtype(r[1])=="K3" and int(c2)==24) else dtype(r[1])
    dpez = [[int(v) for v in L(e)] for e in L(P['DPEZ'][T])]
    ddP  = [e[0] for e in dpez if len(e)>2 and e[2]==1]
    K3   = [i for i,t in types.items() if t=="K3"]
    W    = [i for i,t in types.items() if t=="Wilson"]
    if not (K3 and len(ddP)>=2 and W): return None
    coeur = [K3[0], ddP[0], ddP[1], W[0]]
    for extra in range(1,12):
        if extra in coeur: continue
        basis = coeur+[extra]
        kap = np.zeros((5,5,5))
        for key,val in dict(P['INTNUMS'][T]).items():
            t = tuple(int(x) for x in (key if isinstance(key,tuple)
                       else str(key).strip("()").replace(" ","").split(",")))
            if all(x in basis for x in t):
                a = tuple(basis.index(x) for x in t)
                for p in set(itertools.permutations(a)): kap[p] = float(val)
        if np.count_nonzero(kap) > 3: return kap, basis
    return None

def p_of(kap, t):
    V = np.einsum('ijk,i,j,k',kap,t,t,t)/6.
    if V <= 0: return None
    M = np.einsum('ijk,k',kap,t)
    if abs(np.linalg.det(M)) < 1e-30: return None
    G = -np.linalg.inv(M)/V + np.outer(t,t)/(2*V**2); G = 0.5*(G+G.T)
    tau = np.einsum('ijk,j,k',kap,t,t)/2.
    ev = np.linalg.eigvalsh(G)
    if ev.min() <= 0 or np.any(tau <= 0): return None          # <- LE FILTRE
    if ev.max()/ev.min() > 1e12: return None                   # <- conditionnement
    S = np.linalg.cholesky(G)
    A = np.array([S.T@f for f in [t/(2*V), np.eye(5)[1], np.eye(5)[3]]])
    Q,_ = np.linalg.qr(A.T)
    u,s,_ = np.linalg.svd(np.eye(5)-Q@Q.T)
    B = [np.linalg.solve(S.T, u[:,j]) for j in range(5) if s[j] > 1e-8][:2]
    if len(B) < 2: return None
    B = [b/np.sqrt(b@G@b) for b in B]
    cov = np.eye(5)[0]/tau[0]
    nr = np.sqrt(cov@np.linalg.inv(G)@cov)
    p = float(np.hypot(*[W_GEOM*(cov@b)/nr for b in B]))
    if not np.isfinite(p) or p > W_GEOM*1.0001: return None    # <- LE PLAFOND
    return p

print(f"{'POLY':>6} {'tri':>4} {'h12':>4} {'base':>18} {'adm':>5} {'p_min':>7} "
      f"{'p_max':>7} {'p_med':>7} {'fenetre':>8}")
for pid, T in CIBLES:
    fp = glob.glob(os.path.join(FOLDER, f"POLY_{pid}_*.p"))
    if not fp: print(f"  {pid:>6} : fichier absent"); continue
    with gzip.open(fp[0],"rb") as f: P = pickle.load(f, encoding="latin1")
    r = build(P, T)
    if r is None: print(f"  {pid:6d} {T:4d} : structure absente sur cette triangulation"); continue
    kap, basis = r
    rng = np.random.default_rng(0); v = []
    for _ in range(20000):
        t = np.exp(rng.uniform(-1,4,5))
        Vc = np.einsum('ijk,i,j,k',kap,t,t,t)/6.
        if Vc <= 0: continue
        t *= (1e6/Vc)**(1/3.)
        x = p_of(kap, t)
        if x is not None: v.append(x)
    if len(v) < 20: print(f"  {pid:6d} {T:4d} : {len(v)} points admissibles — trop peu"); continue
    v = np.array(v); fr = float(np.mean((v>=FEN[0])&(v<=FEN[1])))
    print(f"  {pid:6d} {T:4d} {int(P['h12']):4d} {str(basis):>18} {len(v):5d} "
          f"{v.min():7.4f} {v.max():7.4f} {np.median(v):7.4f} {fr:7.1%}")
print("\nCopie-colle ce tableau dans le chat — c'est tout ce dont j'ai besoin.")
