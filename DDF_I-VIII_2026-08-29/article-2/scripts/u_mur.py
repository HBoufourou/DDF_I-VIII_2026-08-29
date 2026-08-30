"""u_mur DERIVE — le vecteur de couplage du mur de domaine. 07/08/2026.
Chaine : T ~ m_chi^4 = M_S^4 = (m_3/2 MPl)^2  =>  T ~ m_3/2^2 (unites de Planck)
et m_3/2 ~ 1/V  =>  ln T = -2 ln V . Il ne reste qu'a exprimer ln V dans notre base."""
import sympy as sp
sq=sp.sqrt
a1=1/(4*sq(7)); b1=-sq(7)/4; b2=-1/sq(35)
print("="*74); print("1. ln(Volume interne, unites de corde) DANS NOTRE BASE")
print("="*74)
print("  espace interne = X5 x intervalle :")
print("    Vol(X5) ~ exp(5(a1 rho + b2 chi))   ;   L ~ exp(b1 rho)")
print("    unites de corde : V = Vol_E * exp(-3 phi/2)   (l_s(Einstein) = e^{phi/4} l_s)")
lnV_phi=sp.Rational(-3,2); lnV_chi=5*b2; lnV_rho=5*a1+b1
print(f"    => ln V = ({lnV_phi}) phi + ({sp.nsimplify(lnV_chi)}) chi + ({sp.nsimplify(lnV_rho)}) rho")
print(f"       numeriquement : ({float(lnV_phi):.5f}, {float(lnV_chi):.5f}, {float(lnV_rho):.5f})")
print("\n"+"="*74); print("2. u_mur = -2 x (vecteur de ln V), en base CANONIQUE (x sqrt2)")
print("="*74)
u=[sp.simplify(-2*x*sq(2)) for x in (lnV_phi,lnV_chi,lnV_rho)]
w=[5*sq(2)/4, sq(sp.Rational(10,7)), 9/(2*sq(14))]
print(f"  u_mur = ({float(u[0]):.4f}, {float(u[1]):.4f}, {float(u[2]):.4f})")
nu=sp.simplify(sum(x**2 for x in u)); nw=sp.simplify(sum(x**2 for x in w))
dot=sp.simplify(sum(a*b for a,b in zip(w,u)))
print(f"\n  >>> |u_mur|^2 = {nu}      <-- TROISIEME THEOREME DE NORME (exact)")
print(f"  >>> |w|^2     = {nw}")
print(f"  >>> w . u     = {dot}       <-- ENTIER")
cos=sp.simplify(dot/(sq(nw)*sq(nu)))
print(f"  >>> cos(angle) = {dot}/sqrt({nw*nu}) = {sp.nsimplify(cos)} = {float(cos):.6f}  ->  {float(sp.deg(sp.acos(cos))):.2f} deg")
print("\n"+"="*74); print("3. c_b : LE MUR EST LA SEULE SOURCE => PLUS AUCUN PARAMETRE LIBRE")
print("="*74)
print("  le puits est creuse par le mur seul (la config symetrique annule les tensions")
print("  de cordes, cf. tadpoles) => M^2 est de rang un LE LONG DE u_mur")
print("  => Phi est l'unique etat massif, n = u_mur/|u_mur|")
cb=sp.simplify(dot/sq(nu))
print(f"\n  >>>>  c_b = (w.u)/|u| = {dot}/sqrt({nu}) = {sp.nsimplify(cb)} = {float(cb):.4f}   [EXACT, aucun parametre]")
print(f"\n  variante (si la pente du puits porte AUSSI le poids du couplage, w+u) :")
s=[sp.simplify(a+b) for a,b in zip(w,u)]
ns=sp.simplify(sum(x**2 for x in s)); cb2=sp.simplify(sum(a*b for a,b in zip(w,s))/sq(ns))
print(f"     |w+u|^2 = {ns} ; c_b = {sp.nsimplify(cb2)} = {float(cb2):.4f}")
print("\n"+"="*74); print("4. VERDICT ET LE MOTIF D'ENSEMBLE")
print("="*74)
print(f"  fenetre du corpus : [0,30 ; 0,61]     obtenu : {float(cb):.3f} (ou {float(cb2):.3f})")
print(f"  => facteur {float(cb)/0.61:.1f} au-dessus du bord haut. EXCLU.")
print("\n  TOUS les c_b derives a ce jour :")
for nom,val in (("LVS fibre (fibre pure)",sp.sqrt(sp.Rational(5,6))),("LVS vanille (volume)",7*sq(6)/12),
                ("mur de domaine",cb),("variante w+u",cb2),("source alignee sur w",sq(6))):
    print(f"    {nom:28s} c_b = {float(val):.4f}")
print("  >>> AUCUN n'est sous 0,91. La fenetre [0,30;0,61] est SOUS TOUT ce que")
print("      le plongement produit — d'un facteur 1,5 a 4. C'est le vrai resultat.")
