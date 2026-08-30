# ══════════════════════════════════════════════════════════════════════════
#  A4 — CALCUL DE |p| SUR TOUS LES CANDIDATS  ·  CELLULE COLAB
#  Lit a4_final.json (déjà chez toi) et calcule la projection pour chacun.
#  Aucun besoin du cône explicite : l'admissibilité est testée par la
#  POSITIVITÉ de la métrique de Kähler, qui est équivalente en pratique.
# ══════════════════════════════════════════════════════════════════════════
import json, itertools, numpy as np
from scipy.optimize import fsolve

HITS = json.load(open("a4_final.json"))
W_GEOM = 1.6956          # |w_geom| du corpus
FEN    = (0.30, 0.61)    # la fenêtre
print(f"{len(HITS)} triangulations à traiter\n")

def kappa_on(basis, intnums, n=5):
    """restreint les INTNUMS aux 5 indices de la base"""
    k = np.zeros((n,n,n)); idx = list(basis)
    for key, val in intnums.items():
        t = tuple(int(x) for x in key.strip("()").replace(" ","").split(","))
        if all(x in idx for x in t):
            a = tuple(idx.index(x) for x in t)
            for p in set(itertools.permutations(a)): k[p] = float(val)
    return k

def metric(kap, t):
    V = np.einsum('ijk,i,j,k',kap,t,t,t)/6.
    M = np.einsum('ijk,k',kap,t)
    if V <= 0 or abs(np.linalg.det(M)) < 1e-30: return None, None, None
    G = -np.linalg.inv(M)/V + np.outer(t,t)/(2*V**2)
    return 0.5*(G+G.T), V, M

def flat_and_p(kap, t, i_small, i_W, i_K3):
    """directions plates (hors volume, τ_s, τ_W) puis |p| pour l'enroulement K3"""
    G, V, M = metric(kap, t)
    if G is None or np.linalg.eigvalsh(G).min() <= 0: return None
    tau = np.einsum('ijk,j,k',kap,t,t)/2.
    if np.any(tau <= 0): return None
    w_,U = np.linalg.eigh(G); S = U@np.diag(np.sqrt(w_))@U.T; Si = U@np.diag(1/np.sqrt(w_))@U.T
    fixes = [t/(2*V), np.eye(5)[i_small], np.eye(5)[i_W]]
    A = np.array([S@f for f in fixes]); Q,_ = np.linalg.qr(A.T)
    u,s,_ = np.linalg.svd(np.eye(5)-Q@Q.T)
    B = [Si@u[:,j] for j in range(5) if s[j] > 1e-8][:2]
    if len(B) < 2: return None
    B = [b/np.sqrt(b@G@b) for b in B]
    cov = np.eye(5)[i_K3]/tau[i_K3]; nr = np.sqrt(cov@np.linalg.inv(G)@cov)
    p = [W_GEOM*(cov@b)/nr for b in B]
    return float(np.hypot(*p))

def scan_point(kap, i_small, i_W, i_K3, Vt=1e6):
    """balaie des points admissibles et renvoie l'intervalle de |p| accessible"""
    rng = np.random.default_rng(0); vals = []
    for _ in range(4000):
        t = np.exp(rng.uniform(-1, 4, 5))
        t *= (Vt/max(np.einsum('ijk,i,j,k',kap,t,t,t)/6., 1e-30))**(1/3.)
        r = flat_and_p(kap, t, i_small, i_W, i_K3)
        if r is not None and np.isfinite(r): vals.append(r)
    return vals

out = []
seen = set()
for h in HITS:
    key = (h["POLYID"], tuple(sorted(h["K3"])), tuple(e[0] for e in h["ddP"]))
    if key in seen: continue            # une géométrie par signature
    seen.add(key)
    K3s = h["K3"]; ddP = [e[0] for e in h["ddP"]]; Ws = h["Wilson"]
    # base physique : K3 + les 2 dP diagonaux + 1 Wilson + un 5e complétant
    coeur = [K3s[0], ddP[0], ddP[1], Ws[0]]
    best = None
    for extra in range(1, 10):
        if extra in coeur: continue
        basis = coeur + [extra]
        kap = kappa_on(basis, h["intnums"])
        if np.count_nonzero(kap) == 0: continue
        i_small, i_W, i_K3 = 1, 3, 0     # positions dans `basis`
        vals = scan_point(kap, i_small, i_W, i_K3)
        if len(vals) > 30:
            best = (basis, np.array(vals)); break
    if best is None:
        out.append({**{k:h[k] for k in ("POLYID","tri","h12","chi")}, "statut":"pas de point admissible"})
        print(f"  POLY {h['POLYID']:5d} tri {h['tri']:3d} : aucun point admissible"); continue
    basis, v = best
    frac = float(np.mean((v >= FEN[0]) & (v <= FEN[1])))
    rec = {"POLYID":h["POLYID"], "tri":h["tri"], "h12":h["h12"], "chi":h["chi"],
           "basis":basis, "n_points":len(v), "p_min":float(v.min()),
           "p_max":float(v.max()), "p_median":float(np.median(v)),
           "fraction_fenetre":frac,
           "fenetre_atteignable": bool(v.min() <= FEN[1] and v.max() >= FEN[0])}
    out.append(rec)
    flag = "  ★ FENÊTRE ATTEIGNABLE" if rec["fenetre_atteignable"] else ""
    print(f"  POLY {h['POLYID']:5d} tri {h['tri']:3d} h12={h['h12']:3d} : "
          f"|p| ∈ [{v.min():.3f} ; {v.max():.3f}] médiane {np.median(v):.3f} "
          f"({frac:.0%} dans la fenêtre){flag}")

print(f"\n{len(out)} géométries distinctes traitées")
ok = [r for r in out if r.get("fenetre_atteignable")]
print(f"{len(ok)} atteignent la fenêtre [0.30 ; 0.61]")
for r in sorted(ok, key=lambda z: -z["fraction_fenetre"])[:15]:
    print(f"   POLY {r['POLYID']} h12={r['h12']} : {r['fraction_fenetre']:.0%} des points, "
          f"médiane {r['p_median']:.3f}")
json.dump(out, open("a4_projections.json","w"), indent=1)
print("\n-> a4_projections.json")
try:
    from google.colab import files; files.download("a4_projections.json")
except Exception: pass
