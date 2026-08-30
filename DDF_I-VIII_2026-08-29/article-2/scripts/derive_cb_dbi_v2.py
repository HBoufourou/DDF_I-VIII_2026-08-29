"""E3a-2 — la matrice cinetique a 3 champs (dilaton, volume X5, radion) et le
vecteur de couplage a 3 composantes. Methode : reduction SEQUENTIELLE
10D --(intervalle)--> 9D --(X5)--> 4D, frame d'Einstein a chaque etape =>
champs canoniques SANS mixage cinetique par construction. 05/08/2026."""
import sympy as sp
sq=sp.sqrt
print("="*74); print("ETAPE 1 — REDUCTION SEQUENTIELLE : les constantes exactes"); print("="*74)
# 10D -> 9D sur l'intervalle (d=9, n=1, D=10) : alpha1^2 = n/(2(d-2)(D-2))
a1=sq(sp.Rational(1,112)); b1=-7*a1
# 9D -> 4D sur X5 (d=4, n=5, D=9)
a2=sq(sp.Rational(5,28)); b2=-2*a2/5
print(f"  intervalle : alpha1 = 1/(4 sqrt7) = {sp.nsimplify(a1)} ; beta1 = -sqrt7/4")
print(f"  X5         : alpha2 = sqrt(5/28)  = {sp.nsimplify(a2)} ; beta2 = -1/sqrt(35)")
print("  metrique finale : ds10^2 = e^{2a1 rho}[ e^{2a2 chi} g4 + e^{2b2 chi} dX5^2 ] + e^{2b1 rho} dy^2")
print("  (rho, chi) canoniques en 4D par construction ; le dilaton phi s'ajoute canonique.")

print("="*74); print("ETAPE 2 — LES POIDS DE LA TENSION D8 (4D x X5, un point sur y)"); print("="*74)
w_phi=sp.Rational(5,4)                      # frame (p-3)/4, p=8
w_chi=sp.simplify(4*a2+5*b2)                # sqrt(det4) x Vol(X5)
w_rho=sp.simplify(9*a1)                     # les 9 dims du volume d'univers portent a1
print(f"  w_phi = 5/4 ; w_chi = 4a2+5b2 = {sp.nsimplify(w_chi)} = {float(w_chi):.4f} ; w_rho = 9a1 = {sp.nsimplify(w_rho)} = {float(w_rho):.4f}")

print("="*74); print("ETAPE 3 — CANONIQUE (x sqrt2) : LE VECTEUR A 3 COMPOSANTES"); print("="*74)
c_phi=sp.simplify(w_phi*sq(2)); c_chi=sp.simplify(w_chi*sq(2)); c_rho=sp.simplify(w_rho*sq(2))
n2=sp.simplify(c_phi**2+c_chi**2+c_rho**2)
print(f"  c_phi = 5 sqrt2/4 = {float(c_phi):.4f} ; c_chi = sqrt(10/7) = {float(c_chi):.4f} ; c_rho = 9/(2 sqrt14) = {float(c_rho):.4f}")
print(f"  |w|^2 = 25/8 + 10/7 + 81/56 = {sp.nsimplify(n2)}   <-- EXACTEMENT 6 : |w| = sqrt6 = {float(sq(n2)):.4f}")
print("  [identite exacte non triviale — a comprendre (somme des poids de Weyl de la")
print("   tension ?), notee comme curiosite, PAS survendue]")

print("="*74); print("ETAPE 4 — CONTROLE DE COHERENCE CONTRE E3a-1 (volume global v6)"); print("="*74)
# direction ou X5 et l'intervalle scalent A L'IDENTIQUE : a1*rho+b2*chi = b1*rho
chi_over_rho=sp.simplify((b1-a1)/b2)        # = 2 sqrt5
nv=sp.Matrix([0,chi_over_rho,1]); nv=nv/nv.norm()
w=sp.Matrix([c_phi,c_chi,c_rho])
proj=sp.simplify((w.T*nv)[0])
print(f"  direction v6 : chi/rho = {sp.nsimplify(chi_over_rho)} ; projection w.n_v6 = {sp.nsimplify(proj)} = {float(proj):.4f}")
print(f"  E3a-1 donnait c_b^volume = 7 sqrt6/12 = {float(7*sq(6)/12):.4f}   ==> {'IDENTIQUE ✓' if sp.simplify(proj-7*sq(6)/12)==0 else 'ECART !'}")

print("="*74); print("ETAPE 5 — LA CONTRAINTE, RADION PROPREMENT EXCLU"); print("="*74)
# le radion du CDD = l'etat propre ~ rho ; Phi vit dans span(phi, chi)
w2=sp.Matrix([c_phi,c_chi]); n2d=sp.simplify(w2.norm())
lo,hi=sp.Rational(30,100),sp.Rational(61,100)
import math
print(f"  sous-espace de Phi = (phi, chi), radion exclu : |w_(phi,chi)| = sqrt(255/56) = {float(n2d):.4f}")
print(f"  fenetre [0.30, 0.61] => |cos(theta)| in [{float(lo/n2d):.3f}, {float(hi/n2d):.3f}]")
print(f"  => Phi a {math.degrees(math.acos(float(hi/n2d))):.1f}-{math.degrees(math.acos(float(lo/n2d))):.1f} degres de la direction de tension : ROBUSTE vs E3a-1 (74-82)")

print("="*74); print("BONUS — LE YUKAWA DU RADION AU NIVEAU ARBRE"); print("="*74)
alpha_rho=sp.simplify(2*c_rho**2)
print(f"  alpha_rho = 2 c_rho^2 = 81/28 = {float(alpha_rho):.3f}")
print("  => la fenetre radion 5-40 um de l'Article I, cherchee avec alpha = O(1),")
print("     recoit un NOMBRE d'arbre : alpha ~ 2,9 (memes hypotheses A1-A3, [tree, declare]).")
print("     Une future carte d'exclusion peut viser cette force precise.")
