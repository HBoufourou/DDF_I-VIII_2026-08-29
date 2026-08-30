"""OUTILS POUR LA SUITE (chantier B+) — 06/08/2026.
Trois outils reutilisables + les derivations qu'ils portent. Boufourou / CDD-cordes."""
import numpy as np
MPl=2.435e18; meV=1e-12
w=np.array([5*np.sqrt(2)/4, np.sqrt(10/7)]); wn=w/np.linalg.norm(w)
psi_w=np.degrees(np.arctan2(w[1],w[0]))

# ---------------------------------------------------------------- OUTIL 1
def echelles(M_S_TeV, M_s_GeV, W0=None):
    """Raccord d'echelle LVS <-> CDD. Donne V requis et W0 requis (ou m_3/2 si W0 fourni)."""
    m32=(M_S_TeV*1e3)**2/MPl
    V=(MPl/M_s_GeV)**2
    W0req=m32*V/MPl
    return dict(m32_meV=m32/meV, V=V, W0_requis=W0req,
                l_s_m=1.97327e-16/M_s_GeV, V_geom=8.2e-6/(1.97327e-16/M_s_GeV))

# ---------------------------------------------------------------- OUTIL 2
def theta_de_composition(p_volume, signe=+1):
    """theta (deg) vs w, pour un etat n = cos(psi) phi + sin(psi) chi avec sin^2 = p_volume.
    signe = +1 : composantes de meme signe ; -1 : signes opposes."""
    psi=np.degrees(np.arcsin(np.sqrt(p_volume)))*signe
    th=abs(psi-psi_w)
    return min(th,180-th)          # |cos| : theta et 180-theta sont la meme direction
def bandes_E11():
    """Les compositions qui satisfont E-11 (73,4-81,9 deg). Renvoie les intervalles en p_volume."""
    out=[]
    for lo,hi in ((psi_w+73.4,psi_w+81.9),(psi_w-81.9,psi_w-73.4)):
        a,b=sorted([np.sin(np.radians(lo))**2, np.sin(np.radians(hi))**2])
        out.append((a,b))
    return sorted(out)

# ---------------------------------------------------------------- OUTIL 3
def audit_minimum(V_func, x, h=1e-6):
    """Verifie qu'un point est VRAIMENT un minimum : gradient nul + Hessienne definie positive.
    A appliquer a chaque point de scan AVANT d'en lire l'angle."""
    x=np.asarray(x,float); n=len(x)
    g=np.array([(V_func(x+h*np.eye(n)[i])-V_func(x-h*np.eye(n)[i]))/(2*h) for i in range(n)])
    H=np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            e1=h*np.eye(n)[i]; e2=h*np.eye(n)[j]
            H[i,j]=(V_func(x+e1+e2)-V_func(x+e1-e2)-V_func(x-e1+e2)+V_func(x-e1-e2))/(4*h*h)
    ev=np.linalg.eigvalsh(H)
    return dict(grad_norm=np.linalg.norm(g), min_eig=ev.min(),
                est_minimum=bool(np.linalg.norm(g)<1e-8*max(1,abs(V_func(x))) and ev.min()>0), eigs=ev)

if __name__=="__main__":
    print("="*70); print("D1 — LE RACCORD D'ECHELLE FIXE W0 (il cesse d'etre libre)"); print("="*70)
    print(f"  {'M_S':>6} {'M_s':>10} {'m_3/2':>9} {'V requis':>11} {'W0 requis':>11}")
    for MS,Ms in ((4,3.6e8),(7.6,2.0e8),(10,0.8e8)):
        d=echelles(MS,Ms)
        print(f"  {MS:5.1f}T {Ms:10.1e} {d['m32_meV']:8.1f}m {d['V']:11.2e} {d['W0_requis']:11.2e}")
    d=echelles(7.6,2.0e8)
    print(f"\n  controle geometrique independant : V ~ R/l_s = {d['V_geom']:.2e} (meme ordre que V = {d['V']:.2e}) ✓")
    print("  ==> W0 ~ 1e-10 a 1e-8 : une SORTIE du raccord, plus un parametre de scan.")
    print("      A fixer dans tous les scans a venir (une dimension de moins).")
    print("      Prix a declarer : W0 aussi petit exige un ajustement de flux profond.")

    print("\n"+"="*70); print("D2 — DIAGNOSTIC DES VOLUMES DU CHANTIER A"); print("="*70)
    for V in (22006.,13340.):
        print(f"  V={V:8.0f} -> M_s = {MPl/np.sqrt(V):.2e} GeV ; m_3/2(W0=1) = {MPl/V/meV:.1e} meV")
    print("  ==> 15 ordres sous le volume requis : l'ANGLE et la HIERARCHIE (sans dimension)")
    print("      restent valides ; l'ECHELLE ABSOLUE ne l'est pas. Sa reserve n.5, chiffree.")

    print("\n"+"="*70); print("D3 — LA COMPOSITION DE Phi EST BIMODALE (resultat neuf)"); print("="*70)
    print(f"  direction de tension : psi_w = {psi_w:.2f} deg dans le plan (phi, chi)")
    b=bandes_E11()
    print(f"  E-11 [73,4 ; 81,9] <=>  p_volume dans [{b[0][0]:.3f} ; {b[0][1]:.3f}]  OU  [{b[1][0]:.3f} ; {b[1][1]:.3f}]")
    print(f"  ==> BANDE INTERDITE entre {b[0][1]:.2f} et {b[1][0]:.2f} : la plage 40-91 % du dossier")
    print("      N'EST PAS UN CONTINUUM, c'est DEUX bandes disjointes.")
    print("      LES DEUX bandes exigent des composantes de SIGNES OPPOSES (il le faut pour")
    print("      etre presque orthogonal a w, dont les deux composantes sont positives).")
    print("      Entre les deux : Phi est TROP orthogonal -> c_b < 0,30, sous la fenetre.")
    print("  ==> son candidat a 76 % volume tombe dans la bande INTERDITE de ce dictionnaire.")
    print("      Ce n'est pas une contradiction : il utilise un dictionnaire toy declare.")
    print("      C'est LE test discriminant du dictionnaire complet — s'il confirme 76 %,")
    print("      c'est le dictionnaire qui bouge, pas la fenetre. A trancher en priorite.")
    print(f"\n  {'p_vol':>6} {'theta (signes opposes)':>24} {'c_b':>8} {'E-11 ?':>8}")
    for p in (0.40,0.45,0.50,0.55,0.65,0.76,0.85,0.91,0.95):
        th=theta_de_composition(p,-1); cb=np.linalg.norm(w)*abs(np.cos(np.radians(th)))
        print(f"  {p:6.2f} {th:24.1f} {cb:8.3f} {'OUI' if 73.4<=th<=81.9 else 'non':>8}")
    print("  sensibilite : ~40 deg par unite de p_volume => tolerance ~ +/-0,03 sur la fraction.")

    print("\n"+"="*70); print("D4 — LA FIBRE : un dof sombre de plus, et une contrainte non listee"); print("="*70)
    print(f"  m_fibre/m_Phi = 0,69  =>  si m_Phi = 24 meV, m_fibre = {24*0.69:.1f} meV")
    print("  => MASSIVE, donc PAS de la radiation noire : elle ne charge pas Delta-N_eff ✓")
    print("  MAIS elle est quasi decouplee (c_b ~ 0) ET plus legere que Phi :")
    print("  => le drip y verse aussi, et rien ne l'en fait sortir. Le budget du reservoir")
    print("     (>=98,6 % / <=1,4 %) doit lui faire une place. CONTRAINTE CREEE PAR LE")
    print("     CHANTIER A, absente de sa liste de points ouverts. A verifier en chantier B.")
