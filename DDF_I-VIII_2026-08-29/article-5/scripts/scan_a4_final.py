# -*- coding: utf-8 -*-
"""A4 FINAL — SCAN AUTORITAIRE DU CQSV. Plus rien à reconstruire.

    python3 scan_a4_final.py CY_Orientifold_database/Data/Complete/h11_5

DÉCOUVERTE : le pickle CQSV contient DÉJÀ tout ce qu'on cherchait à reconstruire —
    DPEZ    = [[idx, n, flag], ...]   les del Pezzo, avec leur degré n et le flag
    DIVS    = [[idx, [h00,h10,h20,h11], chi_top, ?, chi_h], ...]
    SCHERN  = c2.D_i          (le K3 a c2.D = 24)
    FIBRE   = données de fibration
    INTNUMS = les kappa_ijk, par triangulation
Aucun test de rang, aucune matrice de changement de base, aucun ID à faire coïncider.
"""
import pickle, gzip, glob, os, sys, json, itertools
import numpy as np

def load(fp):
    with gzip.open(fp, 'rb') as f: return pickle.load(f)

def divisor_types(divs_T, schern_T):
    """renvoie {index: type} pour une triangulation"""
    out = {}
    for row, c2 in zip(divs_T, schern_T):
        i = int(row[0]); h = list(row[1]); chi_h = float(row[4])
        if h == [1,0,1,20] and int(c2) == 24: out[i] = "K3"
        elif chi_h == 0.0:                    out[i] = "Wilson"
        elif h[1] == 0 and h[2] == 0:         out[i] = f"dP{h[3]-1}"
        else:                                 out[i] = "autre"
    return out

def screen(P):
    """renvoie la liste des triangulations à structure complète"""
    res = []
    nT = len(P['INTNUMS'])
    for T in range(nT):
        types = divisor_types(P['DIVS'][T], P['SCHERN'][T])
        dpez  = [list(map(int, e)) for e in P['DPEZ'][T]]
        ddP7  = [e for e in dpez if e[1] == 7 and e[2] == 1]   # dP7 DIAGONAUX
        K3    = [i for i, t in types.items() if t == "K3"]
        W     = [i for i, t in types.items() if t == "Wilson"]
        if len(K3) >= 1 and len(ddP7) >= 2 and len(W) >= 1:
            res.append({"tri": T, "K3": K3, "ddP7": ddP7, "Wilson": W,
                        "types": {str(k): v for k, v in types.items()},
                        "intnums": {str(k): float(v) for k, v in P['INTNUMS'][T].items()}})
    return res

def run(folder):
    files = sorted(glob.glob(os.path.join(folder, "POLY_*.p")))
    print(f"{len(files)} polytopes\n")
    hits = []
    for fp in files:
        try: P = load(fp)
        except Exception: continue
        good = screen(P)
        if not good: continue
        rec = {"file": os.path.basename(fp), "POLYID": int(P['POLYID']),
               "KSID": int(P['KSID']), "h11": int(P['h11']), "h12": int(P['h12']),
               "chi": int(P['chi']), "triangulations": good}
        hits.append(rec)
        print(f"  ✓ POLY {rec['POLYID']:5d} KSID {rec['KSID']:5d} "
              f"h12={rec['h12']:3d} chi={rec['chi']:5d} : {len(good)} FRST(s) complètes")
    print(f"\n{len(hits)} polytopes à structure complète (K3 + 2 dP7 diagonaux + Wilson)")
    cible = [h for h in hits if h['h12'] == 81 and h['chi'] == -152]
    print(f"   dont à (h12, chi) = (81, -152) — le 4686 de Shukla : {len(cible)}")
    for h in cible: print(f"     POLY {h['POLYID']} / KSID {h['KSID']}")
    autres = [h for h in hits if not (h['h12'] == 81 and h['chi'] == -152)]
    print(f"   AUTRES géométries (candidats neufs) : {len(autres)}")
    for h in autres[:20]: print(f"     POLY {h['POLYID']} h12={h['h12']} chi={h['chi']}")
    json.dump(hits, open("a4_final.json", "w"), indent=1)
    print("\n-> a4_final.json (contient INTNUMS : je peux calculer |p| directement)")

if __name__ == "__main__":
    if len(sys.argv) > 1: run(sys.argv[1])
    else: print(__doc__)
