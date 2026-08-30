"""AUDIT COMPLET DE LA CHAINE c_b — re-derivation independante, maillon par maillon.
Chaque etape est RE-CALCULEE depuis zero (pas recopiee) et confrontee au resultat annonce."""
import sympy as sp
sq=sp.sqrt
OK=lambda b: "OK" if b else "*** ECART ***"
print("="*72); print("M1. REDUCTION KK : les constantes sont-elles celles du frame d'Einstein ?")
print("="*72)
# ds_D^2 = e^{2 a v} ds_d^2 + e^{2 b v} ds_n^2 ; Einstein d-dim : (d-2)a + n b = 0
# canonique : a^2 = n / (2 (d-2)(D-2))
def kk(d,n):
    D=d+n; a=sq(sp.Rational(n,2*(d-2)*(D-2))); b=-(d-2)*a/n; return a,b
a1,b1=kk(9,1); a2,b2=kk(4,5)
print(f"  etape intervalle (d=9,n=1) : a1={sp.nsimplify(a1)} b1={sp.nsimplify(b1)}  ->  attendu 1/(4sqrt7), -sqrt7/4 : {OK(sp.simplify(a1-1/(4*sq(7)))==0 and sp.simplify(b1+sq(7)/4)==0)}")
print(f"  etape X5 (d=4,n=5)         : a2={sp.nsimplify(a2)} b2={sp.nsimplify(b2)}  ->  attendu sqrt(5/28), -1/sqrt35 : {OK(sp.simplify(a2-sq(sp.Rational(5,28)))==0 and sp.simplify(b2+1/sq(35))==0)}")
print(f"  CONTROLE Einstein etape 1 : (d-2)a+nb = {sp.simplify(7*a1+1*b1)} (doit etre 0)")
print(f"  CONTROLE Einstein etape 2 : (d-2)a+nb = {sp.simplify(2*a2+5*b2)} (doit etre 0)")
print("="*72); print("M2. LE PREFACTEUR DE LA TENSION D8 (comptage explicite des exposants)")
print("="*72)
print("  DBI : e^{-phi} sqrt(-det g^S_ind) ; g^S = e^{phi/2} g^E (10D) ; volume d'univers = 9 dims")
p=sp.Integer(8)
chk=sp.simplify(-1+sp.Rational(9,4))  # -phi + (9/2)(phi/2)/... verifions proprement :
# sqrt(det) sur 9 dims du facteur e^{phi/2} : (e^{phi/2})^{9/2} = e^{9phi/4}
w_phi=sp.simplify(-1+sp.Rational(9,4))
print(f"  e^{{-phi}} x e^{{9phi/4}} = e^{{{w_phi} phi}} ;  (p-3)/4 = {sp.Rational(p-3,4)} : {OK(w_phi==sp.Rational(p-3,4))}")
# ds^2 = e^{2a1 r}[e^{2a2 c} g4 + e^{2b2 c} dX5] + e^{2b1 r} dy^2 ; volume d'univers = 4D x X5
w_rho=sp.simplify(9*a1)       # les 9 dims portent toutes e^{a1 r}
w_chi=sp.simplify(4*a2+5*b2)
print(f"  w_rho = 9 a1 = {sp.nsimplify(w_rho)} = {float(w_rho):.4f}")
print(f"  w_chi = 4 a2 + 5 b2 = {sp.nsimplify(w_chi)} = {float(w_chi):.4f}")
print("="*72); print("M3. NORMALISATION CANONIQUE : le facteur sqrt2 est-il le meme pour les 3 ?")
print("="*72)
print("  action reduite : (M^2/2)[R4 - (1/2)(dX)^2] pour CHAQUE X (dilaton et moduli KK)")
print("  => X_c = (M/sqrt2) X  => X = sqrt2 X_c/M  => c_i = sqrt2 w_i  [meme facteur pour les 3] OK")
c_phi=sp.simplify(w_phi*sq(2)); c_chi=sp.simplify(w_chi*sq(2)); c_rho=sp.simplify(w_rho*sq(2))
for n,v,t in (("c_phi",c_phi,5*sq(2)/4),("c_chi",c_chi,sq(sp.Rational(10,7))),("c_rho",c_rho,sp.Rational(9,2)/sq(14))):
    print(f"  {n} = {sp.nsimplify(v)} = {float(v):.4f}   vs annonce {float(t):.4f} : {OK(sp.simplify(v-t)==0)}")
print("="*72); print("M4. LA COMPARAISON A GAROUSI-MYERS : meme OBJET ou coincidence de conventions ?")
print("="*72)
gm=sp.simplify((p-3)/(2*sq(2)))
print(f"  GM (eq.47/eq.6) : (p-3)/(2 sqrt2) = {sp.nsimplify(gm)} = {float(gm):.4f}")
print(f"  nous            : (p-3)/4 x sqrt2 = {sp.nsimplify(c_phi)} = {float(c_phi):.4f}")
print(f"  IDENTITE ALGEBRIQUE : (p-3)/(2sqrt2) == (p-3)sqrt2/4 ? {OK(sp.simplify(gm-(p-3)*sq(2)/4)==0)}")
print("  => ce n'est PAS une coincidence numerique : c'est la MEME expression en p.")
print("  MAIS les deux sqrt2 ont des origines differentes (notre normalisation canonique 4D")
print("  vs leur convention de dilaton 10D). L'accord est reel SI leur phi est canonique ;")
print("  leur L0 = -T_p kappa (h^a_a + (p-3)/(2sqrt2) phi) met h et phi dans le MEME crochet")
print("  avec le meme kappa => meme normalisation. [coherent, mais c'est une INFERENCE]")
print("="*72); print("M5. LE THEOREME |w|^2=6 : re-derivation independante pour 3 branes")
print("="*72)
def vec(p,wrap):
    nt=6-wrap; d1=10-nt
    A1=sq(sp.Rational(nt,2*(d1-2)*8)); A2=sq(sp.Rational(wrap,2*2*(d1-2))); B2=-2*A2/wrap
    wp=sp.Rational(p-3,4); wr=(4+wrap)*A1; wc=4*A2+wrap*B2
    v=[sp.simplify(x*sq(2)) for x in (wp,wr,wc)]
    return v, sp.simplify(sum(x**2 for x in v))
for p_,wr_ in ((8,5),(6,3),(4,1)):
    v,n2=vec(p_,wr_)
    print(f"  D{p_} (enroule {wr_}) : |w|^2 = {n2} : {OK(n2==6)}")
print("  ATTENTION : la formule w_rho=(4+wrap)A1 suppose que TOUT le volume d'univers")
print("  porte le facteur A1 — vrai si la brane est TRANSVERSE aux nt dimensions reduites")
print("  a la premiere etape. Verifie pour la D8 (transverse a l'intervalle seul).")
print("  Pour D6/D4 la premiere etape reduit nt=3,5 dims : la brane doit y etre transverse.")
print("  => le theoreme vaut pour la CLASSE 'brane transverse au premier bloc reduit'. [precise]")
