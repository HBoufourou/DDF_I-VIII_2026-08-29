"""CELLULE DE SCAN h11=5 — del Pezzo diagonal + budget D-terme + directions plates.
A lancer la ou l'acces aux donnees existe (Colab / HuggingFace / cytools en ligne).
08/08/2026"""
import numpy as np, itertools

# ---------------------------------------------------------------- TESTS
def del_pezzo_diagonal(kappa, n):
    """kappa[i,j,k] : nombres d'intersection. Renvoie la liste des diviseurs diagonaux.
    D_s est diagonal si kappa_sss != 0 et kappa_sij = 0 pour tout (i,j) != (s,s)."""
    out=[]
    for s in range(n):
        if abs(kappa[s,s,s])<1e-9: continue
        ok=True
        for i in range(n):
            for j in range(n):
                if (i,j)!=(s,s) and abs(kappa[s,i,j])>1e-9: ok=False; break
            if not ok: break
        if ok: out.append(s)
    return out

def metrique_log(kappa, tau, n):
    """G_ij en coordonnees log, depuis K = -2 ln V. Ne suppose PAS la diagonalite."""
    # V a partir des 2-cycles serait plus exact ; ici approximation par la forme de volume
    # fournie par l'appelant. Placeholder : a remplacer par le calcul cytools exact.
    raise NotImplementedError("remplacer par cy.compute_kahler_metric() de cytools")

def directions_plates(G, dir_volume, dir_blowups, dir_Dterme=None):
    """Sous-espace orthogonal (au sens de G) au volume, aux blow-ups, et au D-terme.
    Renvoie une base orthonormee-G des directions plates restantes."""
    fixes=[dir_volume]+list(dir_blowups)+([dir_Dterme] if dir_Dterme is not None else [])
    n=G.shape[0]; B=[]
    for e in np.eye(n):
        v=e.copy()
        for f in fixes+B:
            v=v-(v@G@f)/(f@G@f)*f
        if np.sqrt(abs(v@G@v))>1e-8:
            B.append(v/np.sqrt(v@G@v))
    return B

def projections(base_plate, w_geom_dir, norme_w_geom=1.6956):
    """p_i = projection du vecteur de tension sur chaque direction plate."""
    return [norme_w_geom*(w_geom_dir@G@b) for b in base_plate]  # G dans la portee de l'appelant

def cb_range(ps):
    R=np.sqrt(sum(p*p for p in ps)); return 0.0, R

def fraction_fenetre(ps, lo=0.30, hi=0.61, N=200000):
    """Fraction des angles donnant c_b dans la fenetre (2 directions plates)."""
    if len(ps)<2: return 0.0
    t=np.random.default_rng(0).uniform(0,2*np.pi,N)
    cb=np.abs(np.cos(t)*ps[0]+np.sin(t)*ps[1])
    return float(np.mean((cb>=lo)&(cb<=hi)))

# ---------------------------------------------------------------- PIPELINE
def evalue_candidat(cy, dir_Dterme=None):
    """cy : objet cytools. Renvoie un verdict complet."""
    n = cy.h11()
    kappa = cy.intersection_numbers(in_basis=True, format="dense")
    dp = del_pezzo_diagonal(kappa, n)
    if not dp:
        return dict(verdict="REJET : aucun del Pezzo diagonal (runaway garanti)")
    G = cy.compute_kahler_metric()           # a evaluer au point LVS
    dir_vol = np.ones(n)                      # a remplacer par grad(ln V) exact
    plates = directions_plates(G, dir_vol, [np.eye(n)[s] for s in dp], dir_Dterme)
    return dict(n_del_pezzo=len(dp), n_plates=len(plates),
                verdict=("OK : 2 directions plates" if len(plates)>=2
                         else f"REJET : {len(plates)} direction(s) plate(s) apres prelevement"))

# ---------------------------------------------------------------- BOUCLE
if __name__ == "__main__":
    from cytools import fetch_polytopes
    resultats=[]
    for p in fetch_polytopes(h11=5, lattice="N", limit=500):
        cy = p.triangulate().get_cy()
        if abs(cy.chi())>960 or cy.h21()<=cy.h11(): continue
        r = evalue_candidat(cy)
        if r.get("n_plates",0)>=2: resultats.append((cy, r))
    print(len(resultats), "candidats gardent 2 directions plates apres prelevement D-terme")
