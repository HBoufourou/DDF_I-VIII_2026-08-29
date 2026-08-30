# ══════════════════════════════════════════════════════════════════════════
#  A4 — SCAN CQSV  ·  CELLULE COLAB UNIQUE  ·  copier-coller tel quel
#  Ne demande RIEN d'autre. Clone la base si besoin, trouve les dossiers
#  tout seul, tourne en 2 modes : DIAGNOSTIC puis SCAN.
# ══════════════════════════════════════════════════════════════════════════
import os, sys, glob, gzip, pickle, json, subprocess

REPO = "CY_Orientifold_database"
H11  = 5                    # mettre 6 pour scanner h11=6 ensuite
MODE = "diagnostic"         # "diagnostic" d'abord, puis repasser à "scan"

# ---------- 1. la base -----------------------------------------------------
if not os.path.isdir(REPO):
    print("clonage de la base (quelques minutes)...")
    subprocess.run(["git", "clone", "--depth", "1",
                    "https://github.com/AndreasSchachner/CY_Orientifold_database.git"],
                   check=False)

# ---------- 2. trouver le dossier, quelle que soit l'arborescence -----------
cands = [p for p in glob.glob(f"**/h11_{H11}", recursive=True) if os.path.isdir(p)]
cands += [p for p in glob.glob(f"**/h11={H11}", recursive=True) if os.path.isdir(p)]
FOLDER = None
for c in cands:
    if glob.glob(os.path.join(c, "POLY_*.p")):
        FOLDER = c; break
if FOLDER is None:
    print("!! dossier introuvable. Dossiers vus :")
    for c in sorted(set(cands))[:20]: print("   ", c, len(os.listdir(c)), "fichiers")
    print("\nArborescence du dépôt :")
    for r, ds, fs in os.walk(REPO):
        if fs[:1]: print("   ", r, "->", len(fs), "fichiers", fs[:2])
        if r.count(os.sep) > 3: ds[:] = []
    raise SystemExit("Corrige FOLDER à la main puis relance.")
FILES = sorted(glob.glob(os.path.join(FOLDER, "POLY_*.p")))
print(f"dossier : {FOLDER}\nfichiers : {len(FILES)}\n")

# ---------- 3. lecture robuste --------------------------------------------
def load(fp):
    with gzip.open(fp, "rb") as f:
        return pickle.load(f, encoding="latin1")

def as_list(x):
    """DIVS/DPEZ arrivent parfois en ndarray, parfois en liste"""
    try:    return x.tolist()
    except Exception: return list(x)

def classify(divs_T, schern_T):
    """{index: type} — K3 / Wilson / dP_n, à partir des Hodge + c2"""
    out = {}
    for row, c2 in zip(as_list(divs_T), as_list(schern_T)):
        r = as_list(row)
        i = int(r[0]); h = [int(v) for v in as_list(r[1])]; chi_h = float(r[4])
        if h[:4] == [1,0,1,20] and int(c2) == 24: out[i] = "K3"
        elif chi_h == 0.0:                        out[i] = "Wilson"
        elif h[1] == 0 and h[2] == 0:             out[i] = f"dP{h[3]-1}"
        else:                                      out[i] = "autre"
    return out

# ---------- 4. MODE DIAGNOSTIC : on regarde UN fichier ---------------------
if MODE == "diagnostic":
    P = load(FILES[0])
    print("=== clés du pickle ==="); print(list(P.keys()))
    print(f"\nPOLYID {P['POLYID']}  KSID {P['KSID']}  h11 {P['h11']}  h12 {P['h12']}  chi {P['chi']}")
    print(f"nombre de triangulations : {len(P['INTNUMS'])}")
    print("\n=== DIVS[0] ==="); print(as_list(P['DIVS'][0])[:4])
    print("\n=== SCHERN[0] ==="); print(as_list(P['SCHERN'][0]))
    print("\n=== DPEZ[0] ==="); print(as_list(P['DPEZ'][0]))
    print("\n=== classification obtenue ===")
    for k, v in classify(P['DIVS'][0], P['SCHERN'][0]).items(): print(f"   diviseur {k} -> {v}")
    print("\n>>> Si cette sortie a du sens, remets MODE = 'scan' et relance.")
    raise SystemExit

# ---------- 5. MODE SCAN ---------------------------------------------------
hits = []
for n, fp in enumerate(FILES):
    if n % 200 == 0: print(f"  ... {n}/{len(FILES)}")
    try: P = load(fp)
    except Exception: continue
    for T in range(len(P['INTNUMS'])):
        try:
            types = classify(P['DIVS'][T], P['SCHERN'][T])
            dpez  = [[int(v) for v in as_list(e)] for e in as_list(P['DPEZ'][T])]
        except Exception: continue
        ddP7 = [e for e in dpez if len(e) >= 3 and e[1] == 7 and e[2] == 1]
        K3   = [i for i, t in types.items() if t == "K3"]
        W    = [i for i, t in types.items() if t == "Wilson"]
        if len(K3) >= 1 and len(ddP7) >= 2 and len(W) >= 1:
            hits.append({
                "file": os.path.basename(fp), "POLYID": int(P['POLYID']),
                "KSID": int(P['KSID']), "h11": int(P['h11']), "h12": int(P['h12']),
                "chi": int(P['chi']), "tri": T, "K3": K3, "ddP7": ddP7, "Wilson": W,
                "types": {str(k): v for k, v in types.items()},
                "intnums": {str(k): float(v) for k, v in dict(P['INTNUMS'][T]).items()},
            })
            print(f"  ✓ POLY {P['POLYID']} tri {T} : h12={P['h12']} chi={P['chi']} "
                  f"K3={K3} ddP7={[e[0] for e in ddP7]} W={W}")

print(f"\n{len(hits)} triangulations à structure complète")
cible  = [h for h in hits if h['h12'] == 81 and h['chi'] == -152]
autres = [h for h in hits if not (h['h12'] == 81 and h['chi'] == -152)]
print(f"  à (h12,chi) = (81,-152)  [le 4686 de Shukla] : {len(cible)}"
      f"  -> POLYIDs {sorted(set(h['POLYID'] for h in cible))}")
print(f"  AUTRES géométries [candidats neufs]          : {len(autres)}"
      f"  -> POLYIDs {sorted(set(h['POLYID'] for h in autres))}")

with open("a4_final.json", "w") as f: json.dump(hits, f, indent=1)
print("\n-> a4_final.json écrit")
try:
    from google.colab import files; files.download("a4_final.json")
except Exception: pass
