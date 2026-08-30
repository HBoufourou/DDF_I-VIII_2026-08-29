# ══════════════════════════════════════════════════════════════════════════
#  A4 — SCAN CQSV  ·  CELLULE COLAB  ·  copier-coller, rien à modifier
#  Corrigé : pas de SystemExit (plus de traceback), et le critère porte sur
#  les del Pezzo DIAGONAUX de degré quelconque — pas seulement dP7.
# ══════════════════════════════════════════════════════════════════════════
import os, glob, gzip, pickle, json, warnings, subprocess
from collections import Counter
warnings.filterwarnings("ignore")

H11 = 5          # passer à 6 ensuite si tu veux élargir
REPO = "CY_Orientifold_database"

if not os.path.isdir(REPO):
    subprocess.run(["git","clone","--depth","1",
        "https://github.com/AndreasSchachner/CY_Orientifold_database.git"], check=False)

FOLDER = None
for pat in (f"**/h11_{H11}", f"**/h11={H11}"):
    for c in glob.glob(pat, recursive=True):
        if os.path.isdir(c) and glob.glob(os.path.join(c,"POLY_*.p")): FOLDER = c; break
    if FOLDER: break
FILES = sorted(glob.glob(os.path.join(FOLDER,"POLY_*.p"))) if FOLDER else []
print(f"dossier : {FOLDER}\nfichiers : {len(FILES)}\n")

def load(fp):
    with gzip.open(fp,"rb") as f: return pickle.load(f, encoding="latin1")
def L(x):
    try: return x.tolist()
    except Exception: return list(x)

def classify(divs_T, schern_T):
    out = {}
    for row, c2 in zip(L(divs_T), L(schern_T)):
        r = L(row); i = int(r[0]); h = [int(v) for v in L(r[1])]; chi_h = float(r[4])
        if h[:4] == [1,0,1,20] and int(c2) == 24: out[i] = "K3"
        elif chi_h == 0.0:                        out[i] = "Wilson"
        elif h[1] == 0 and h[2] == 0:             out[i] = f"dP{h[3]-1}"
        else:                                      out[i] = "autre"
    return out

hits, flags, nT_tot = [], Counter(), 0
for n, fp in enumerate(FILES):
    if n % 400 == 0: print(f"  ... {n}/{len(FILES)}")
    try: P = load(fp)
    except Exception: continue
    for T in range(len(P["INTNUMS"])):
        nT_tot += 1
        try:
            types = classify(P["DIVS"][T], P["SCHERN"][T])
            dpez  = [[int(v) for v in L(e)] for e in L(P["DPEZ"][T])]
        except Exception: continue
        for e in dpez: flags[e[2] if len(e) > 2 else None] += 1
        # ---- critère : >= 2 del Pezzo DIAGONAUX (flag = 1), degré quelconque
        ddP = [e for e in dpez if len(e) > 2 and e[2] == 1]
        K3  = [i for i,t in types.items() if t == "K3"]
        W   = [i for i,t in types.items() if t == "Wilson"]
        if len(K3) >= 1 and len(ddP) >= 2 and len(W) >= 1:
            hits.append({"file":os.path.basename(fp), "POLYID":int(P["POLYID"]),
                "KSID":int(P["KSID"]), "h11":int(P["h11"]), "h12":int(P["h12"]),
                "chi":int(P["chi"]), "tri":T, "K3":K3, "Wilson":W,
                "ddP":ddP, "degres_ddP":[e[1] for e in ddP],
                "types":{str(k):v for k,v in types.items()},
                "intnums":{str(k):float(v) for k,v in dict(P["INTNUMS"][T]).items()}})
            print(f"  ✓ POLY {P['POLYID']} tri {T} : h12={P['h12']} chi={P['chi']} "
                  f"K3={K3} ddP={[(e[0],'dP%d'%e[1]) for e in ddP]} W={W}")

print(f"\n{nT_tot} triangulations parcourues")
print(f"répartition des flags DPEZ : {dict(flags)}   (1 = diagonal)")
print(f"{len(hits)} triangulations à structure complète")
if hits:
    cible  = [h for h in hits if h["h12"]==81 and h["chi"]==-152]
    autres = [h for h in hits if not (h["h12"]==81 and h["chi"]==-152)]
    print(f"  (h12,chi)=(81,-152) [le 4686 de Shukla] : {len(cible)} "
          f"-> POLYIDs {sorted(set(h['POLYID'] for h in cible))}")
    print(f"  AUTRES géométries  [candidats neufs]    : {len(autres)} "
          f"-> POLYIDs {sorted(set(h['POLYID'] for h in autres))}")
    print(f"  degrés de dP diagonaux rencontrés : "
          f"{Counter(d for h in hits for d in h['degres_ddP'])}")
    with open("a4_final.json","w") as f: json.dump(hits, f, indent=1)
    print("\n-> a4_final.json écrit")
    try:
        from google.colab import files; files.download("a4_final.json")
    except Exception: pass
else:
    print("\nAucun résultat. Colle-moi la répartition des flags ci-dessus :")
    print("si tous les flags valent 0, le flag n'est pas la diagonalité et on change de critère.")
