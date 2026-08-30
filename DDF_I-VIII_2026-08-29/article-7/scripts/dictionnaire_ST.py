"""POINT 1 — LE DICTIONNAIRE (S,T) <-> (phi, chi, rho). 05/08/2026.
Principe : les DEUX bases sont canoniques, donc le passage est une ISOMETRIE.
Il suffit d'identifier les DIRECTIONS geometriques de S et T dans la base
canonique (phi, chi, rho) de la reduction sequentielle, puis d'orthonormaliser."""
import numpy as np
a1=1/(4*np.sqrt(7)); b1=-np.sqrt(7)/4; a2=np.sqrt(5/28); b2=-1/np.sqrt(35)
print(f"constantes de la reduction : a1={a1:.4f} b1={b1:.4f} a2={a2:.4f} b2={b2:.4f}\n")
print("=== 1. LES DEUX QUANTITES GEOMETRIQUES QUI DEFINISSENT S ET T ===")
print("  Re S  <->  couplage de jauge de la D8 : 1/g^2 ~ e^{-phi_10} Vol(X5)")
print("  Re T  <->  volume interne total : V6 = Vol(X5) x L_intervalle")
# phi_10 = sqrt2 phi_c ; ln Vol(X5) = 5(a1 rho + b2 chi) ; ln L = b1 rho
gS=np.array([-np.sqrt(2), 5*b2, 5*a1])        # grad ln(1/g^2) dans (phi, chi, rho)
gT=np.array([0.0,        5*b2, 5*a1+b1])      # grad ln V6
print(f"  grad ln(1/g^2) = ({gS[0]:+.3f}, {gS[1]:+.3f}, {gS[2]:+.3f})")
print(f"  grad ln V6     = ({gT[0]:+.3f}, {gT[1]:+.3f}, {gT[2]:+.3f})")
print("\n=== 2. ORTHONORMALISATION (Gram-Schmidt) : les directions S et T ===")
nS=gS/np.linalg.norm(gS)
t=gT-(gT@nS)*nS; nT=t/np.linalg.norm(t)
print(f"  n_S = ({nS[0]:+.4f}, {nS[1]:+.4f}, {nS[2]:+.4f})   [phi, chi, rho]")
print(f"  n_T = ({nT[0]:+.4f}, {nT[1]:+.4f}, {nT[2]:+.4f})")
print(f"  verif orthogonalite n_S.n_T = {nS@nT:+.2e}")
print("\n=== 3. LE DICTIONNAIRE APPLIQUE AU VIDE RACETRACK ===")
nl_sigma,nl_tau=0.616,-0.788                   # etat leger trouve en (sigma, tau)
n_light=nl_sigma*nS+nl_tau*nT
print(f"  etat leger du racetrack (0,616 sigma - 0,788 tau) devient, dans (phi,chi,rho) :")
print(f"    n = ({n_light[0]:+.4f}, {n_light[1]:+.4f}, {n_light[2]:+.4f})")
w=np.array([5*np.sqrt(2)/4, np.sqrt(10/7), 9/(2*np.sqrt(14))]); wh=w/np.linalg.norm(w)
th=np.degrees(np.arccos(min(1,abs(n_light@wh))))
cb=abs(n_light@w)
print(f"\n=== 4. LE RESULTAT : L'ANGLE ET c_b DE CE VIDE ===")
print(f"  angle avec la direction de tension : theta = {th:.1f} deg")
print(f"  fenetre CDD exigee : 73,4 - 81,9 deg  ->  {'DANS LA FENETRE' if 73.4<=th<=81.9 else 'HORS FENETRE'}")
print(f"  c_b = |n . w| = {cb:.4f}   (fenetre corpus [0,30 ; 0,61])")
print(f"  contenu radion de l'etat : {n_light[2]**2:.1%}  (doit rester faible : le radion est un AUTRE champ)")
print("\n=== 5. CE QUE CE PREMIER VIDE DIT ===")
if not (73.4<=th<=81.9):
    print("  ce vide-ci N'EST PAS un vide du CDD : il echoue E-11 (comme sa hierarchie 0,71")
    print("  l'annonçait deja). C'est le resultat ATTENDU d'un premier point de parametres")
    print("  non scanne — et la chaine complete (K,W -> minimum -> Hessienne -> dictionnaire")
    print("  -> theta -> c_b) FONCTIONNE maintenant de bout en bout. Le scan peut commencer.")
print("\n  CONVENTION DECLAREE : Re S ~ e^{-phi}Vol(X5) (couplage D8) et Re T ~ V6.")
print("  Une autre affectation (T = cycle particulier) change n_T et donc theta :")
print("  c'est LE point ou une construction explicite doit fixer le choix. [Open, declare]")
