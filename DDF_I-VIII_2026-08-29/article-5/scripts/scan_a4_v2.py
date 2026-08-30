# -*- coding: utf-8 -*-
"""A4 v2 — DEUX CORRECTIFS. À relancer sur le bucket NID.
    python3 scan_a4_v2.py NID_h11_5/NID/h11=5

CORRECTIF 1 : lire divisors_hodge sur TOUS les diviseurs toriques (0..N-1),
              pas seulement sur les 5 de la base. Shukla décrit 4686 par SIX
              diviseurs alors que h11 = 5 : au moins un est hors base.
CORRECTIF 2 : chercher la diagonalité sur les CLASSES v du réseau, pas
              seulement sur v = e_i. rang(κ(v,·,·)) = 1 est invariant, mais
              il faut balayer v.

En sortie, chaque FRST porte aussi le tenseur κ des candidats retenus, pour
que le calcul de |p| puisse suivre sans relire le bucket.
"""
import os, re, sys, json, ast, itertools
import numpy as np

def P(v):
    if not isinstance(v, str): return v
    try: return ast.literal_eval(v)
    except Exception: return v

def hodge(kk):
    m = re.search(r'H:(\d+),(\d+)', str(kk))
    return (int(m.group(1)), int(m.group(2))) if m else None

def dtype(h):
    h00, h10, h20, h11 = h[:4]
    if h10 == 0 and h20 == 1 and h11 == 20: return "K3"
    if h10 == 1 and h20 == 0:               return "Wilson"
    if h10 == 0 and h20 == 0:               return f"dP{h11-1}"
    return f"other{tuple(h[:4])}"

def kap_full(triple, n):
    """tenseur dense sur TOUS les indices toriques"""
    k = np.zeros((n,n,n))
    for key, val in triple.items():
        t = tuple(key) if not isinstance(key, str) else ast.literal_eval(key)
        for p in set(itertools.permutations(t)): k[p] = float(val)
    return k

def kap_basis(triple, idx):
    n = len(idx); k = np.zeros((n,n,n))
    for key, val in triple.items():
        t = tuple(key) if not isinstance(key, str) else ast.literal_eval(key)
        if all(x in idx for x in t):
            a = tuple(idx.index(x) for x in t)
            for p in set(itertools.permutations(a)): k[p] = float(val)
    return k

def diagonal_classes(kap, n, rng=2, tol=1e-9):
    """CORRECTIF 2 : toutes les classes primitives v telles que rang(κ(v,·,·)) = 1"""
    out = []
    for c in itertools.product(range(-rng, rng+1), repeat=n):
        if not any(c): continue
        if np.gcd.reduce([abs(x) for x in c]) != 1: continue
        v = np.array(c, float)
        M = np.einsum('ijk,i->jk', kap, v)
        if np.linalg.matrix_rank(M, tol=tol) != 1: continue
        if abs(np.einsum('jk,j,k', M, v, v)) > tol: out.append(list(map(int, c)))
    return out

def screen(f):
    hh = hodge(f.get("KK_data", ""))
    if not hh: return None
    h21, h11 = hh
    dinq = P(f.get("divisor_indenp")); dh = P(f.get("divisors_hodge")); tri = P(f.get("triple_inters"))
    if not (dinq and dh and tri): return None
    idx = [int(re.sub(r'\D', '', str(z))) for z in dinq]
    dh  = [P(x) for x in dh]
    N   = len(dh)
    # --- CORRECTIF 1 : tous les diviseurs toriques
    types_all  = [dtype(h) for h in dh]
    types_base = [types_all[i] for i in idx]
    # --- CORRECTIF 2 : classes diagonales dans la base
    kb   = kap_basis(tri, idx)
    diag = diagonal_classes(kb, len(idx))
    dP_base = [i for i,t in enumerate(types_base) if t.startswith("dP")]
    return {
        "h11": h11, "h21": h21, "chi": 2*(h11-h21), "basis": idx,
        "types_base": types_base, "types_ALL": types_all,
        "nK3_base": types_base.count("K3"), "nK3_ALL": types_all.count("K3"),
        "nW_ALL": types_all.count("Wilson"),
        "diag_classes": diag, "n_diag": len(diag),
        "dP_in_base": [types_base[i] for i in dP_base],
        "kappa_basis": kb.tolist(),
        "COMPLET": (types_all.count("K3") >= 1 and len(diag) >= 2
                    and types_all.count("Wilson") >= 1),
    }

def run(folder):
    hits = []
    files = sorted(x for x in os.listdir(folder) if x.endswith('.json'))
    for fn in files:
        try: data = json.load(open(os.path.join(folder, fn)))
        except Exception: continue
        for t, frst in enumerate(data if isinstance(data, list) else [data]):
            r = screen(frst)
            if not r: continue
            r["file"], r["tri"] = fn, t
            if r["COMPLET"]:
                hits.append(r)
                print(f"  ✓ {fn}#{t} : h21={r['h21']} chi={r['chi']} "
                      f"K3(all)={r['nK3_ALL']} diag={r['n_diag']} W={r['nW_ALL']}")
    print(f"\n{len(hits)} FRSTs à structure complète")
    json.dump(hits, open("a4_v2_complets.json", "w"), indent=1)
    print("-> a4_v2_complets.json (contient kappa_basis pour le calcul de |p|)")

if __name__ == "__main__":
    if len(sys.argv) > 1: run(sys.argv[1])
    else: print(__doc__)
