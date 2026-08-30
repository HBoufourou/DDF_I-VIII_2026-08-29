# -*- coding: utf-8 -*-
"""A4 — SCAN DU BUCKET NID PAR STRUCTURE, PAS PAR IDENTIFIANT.
    python3 scan_a4.py  NID_h11_5/NID/h11=5

Écrit pour les champs RÉELS de ton bucket : triple_inters, divisor_indenp,
divisors_hodge, KK_data. Le bucket REF n'est pas nécessaire.

POURQUOI PAR STRUCTURE : le polyid 4686 du NID a H:121,5 (h21=121, chi=-232),
alors que le 4686 de Shukla a h21=81, chi=-152. Les numérotations diffèrent.
On cherche donc la GÉOMÉTRIE, pas le numéro.
"""
import os, re, sys, json, ast, itertools
import numpy as np

# ----------------------------------------------------------------- PARSING
def parse_field(v):
    """les champs sont souvent des chaînes contenant du Python/Mathematica"""
    if not isinstance(v, str): return v
    try: return ast.literal_eval(v)
    except Exception: return v

def hodge_from_KK(kk):
    """'... H:121,5 [232]' -> (h21, h11, chi)"""
    m = re.search(r'H:(\d+),(\d+)', str(kk))
    if not m: return None
    h21, h11 = int(m.group(1)), int(m.group(2))
    return h21, h11, 2*(h11 - h21)

def divisor_type(hv):
    """[h00,h10,h20,h11] -> étiquette"""
    h00, h10, h20, h11 = hv[:4]
    chi_h = h00 - h10 + h20
    if h20 == 1 and h11 == 20 and h10 == 0:      return "K3"
    if h10 == 1 and h20 == 0:                    return "Wilson"      # chi_h = 0
    if h10 == 0 and h20 == 0:                    return f"dP{h11-1}"  # del Pezzo de rang h11
    return f"other({h00},{h10},{h20},{h11})"

def kappa_in_basis(triple, idx, n):
    """restreint triple_inters aux indices de base -> tenseur dense n^3"""
    kap = np.zeros((n,n,n))
    for key, val in triple.items():
        k = tuple(key) if not isinstance(key, str) else ast.literal_eval(key)
        if all(x in idx for x in k):
            a, b, c = (idx.index(x) for x in k)
            for p in set(itertools.permutations((a,b,c))): kap[p] = float(val)
    return kap

# ----------------------------------------------------------------- LES TESTS
def is_diagonal(kap, i, tol=1e-9):
    """dP diagonal : rang(kappa(e_i,.,.)) == 1 et kappa_iii != 0 — INVARIANT DE BASE"""
    M = kap[i]
    return np.linalg.matrix_rank(M, tol=tol) == 1 and abs(kap[i,i,i]) > tol

def screen_frst(f):
    """renvoie un verdict pour une FRST, ou None si les champs manquent"""
    kk = hodge_from_KK(f.get("KK_data",""))
    if not kk: return None
    h21, h11, chi = kk
    dinq = parse_field(f.get("divisor_indenp"))
    dh   = parse_field(f.get("divisors_hodge"))
    tri  = parse_field(f.get("triple_inters"))
    if not (dinq and dh and tri): return None
    idx = [int(re.sub(r'\D','',str(z))) for z in dinq]
    dh  = [parse_field(x) if isinstance(x,str) else x for x in dh]
    types = [divisor_type(dh[i]) for i in idx]
    kap = kappa_in_basis(tri, idx, len(idx))
    diag = [k for k in range(len(idx)) if is_diagonal(kap, k)]
    dP_diag = [k for k in diag if types[k].startswith("dP")]
    return {"h11":h11, "h21":h21, "chi":chi, "basis":idx, "types":types,
            "dP_diagonaux":[types[k] for k in dP_diag],
            "n_K3":types.count("K3"), "n_Wilson":types.count("Wilson"),
            "CANDIDAT": (h11>=4 and h21>h11 and abs(chi)<=960
                         and types.count("K3")>=1 and len(dP_diag)>=2
                         and types.count("Wilson")>=1),
            "CIBLE_SHUKLA": (h11==5 and h21==81 and chi==-152)}

# ----------------------------------------------------------------- BOUCLE
def run(folder):
    files = sorted(f for f in os.listdir(folder) if f.endswith('.json'))
    print(f"{len(files)} fichiers, chacun contenant plusieurs FRSTs\n")
    hits, cibles, nfrst = [], [], 0
    for fn in files:
        try: data = json.load(open(os.path.join(folder, fn)))
        except Exception: continue
        for t, frst in enumerate(data if isinstance(data, list) else [data]):
            nfrst += 1
            r = screen_frst(frst)
            if not r: continue
            r["file"], r["tri"] = fn, t
            if r["CIBLE_SHUKLA"]:
                cibles.append(r)
                print(f"  ★ CIBLE  {fn}#{t} : h21=81, chi=-152, types={r['types']}")
            if r["CANDIDAT"]:
                hits.append(r)
                print(f"    cand.  {fn}#{t} : chi={r['chi']}, K3x{r['n_K3']}, "
                      f"dP diag {r['dP_diagonaux']}, Wx{r['n_Wilson']}")
    print(f"\n{nfrst} FRSTs scannées")
    print(f"  geometries a (h11,h21,chi) = (5,81,-152)  : {len(cibles)}")
    print(f"  CANDIDATS (structure complete)             : {len(hits)}")
    json.dump({"cibles":cibles, "candidats":hits}, open("a4_resultats.json","w"),
              indent=1, default=str)
    print("-> a4_resultats.json")

if __name__ == "__main__":
    if len(sys.argv) > 1: run(sys.argv[1])
    else: print(__doc__)
