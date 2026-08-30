"""E3-FINAL — le maximum derivable sur c_b et toutes les valeurs numeriques. 05/08/2026.
1) THEOREME |w|^2 = 6 (teste exactement pour 3 configurations de branes)
2) c_b^max exact, alpha_b^max, exclusion des benchmarks hauts du corpus
3) la COMPOSITION de Phi (fraction volume/dilaton) imposee par la fenetre
4) la distribution de c_b sous la mesure de cordes (le max derivable sans stabilisation)"""
import sympy as sp, numpy as np
sq=sp.sqrt
def config(p,wrap):
    """reduction sequentielle : transverse d'abord, puis les dims enroulees."""
    nt=6-wrap                                   # transverses internes
    d1=10-nt
    a1=sq(sp.Rational(nt,2*(d1-2)*8)); b1=-(d1-2)*a1/nt
    a2=sq(sp.Rational(wrap,2*2*(d1-2))); b2=-2*a2/wrap
    w_phi=sp.Rational(p-3,4)
    w_rho=(4+wrap)*a1                            # tout le volume d'univers porte a1
    w_chi=4*a2+wrap*b2
    v=[sp.simplify(x*sq(2)) for x in (w_phi,w_rho,w_chi)]
    return v, sp.simplify(sum(x**2 for x in v))
print("=== 1. LE THEOREME : |w|^2 = 6, INDEPENDANT de p et de l'enroulement ===")
for p,wrap in ((8,5),(6,3),(4,1)):
    v,n2=config(p,wrap)
    print(f"  D{p} enroulant {wrap} dims : c = ({', '.join(str(sp.nsimplify(x)) for x in v)}) ; |w|^2 = {n2}")
print("  => la 'curiosite' d'E3a-2 est un THEOREME de la reduction : la norme du vecteur")
print("     de couplage de la tension vaut sqrt6 pour TOUTE Dp. [R11 A FAIRE : poids")
print("     conforme universel des termes de tension — probablement connu (litterature")
print("     des potentiels exponentiels / swampland) ; a verifier AVANT tout claim.]")
print()
print("=== 2. LE MAXIMUM DE c_b (exact) ET SES CONSEQUENCES ===")
c_phi=5*sq(2)/4; c_chi=sq(sp.Rational(10,7)); c_rho=sp.Rational(9,2)/sq(14)
cbmax=sp.simplify(sq(c_phi**2+c_chi**2))       # radion exclu
w1=5.48
print(f"  c_b^max = |w_(phi,chi)| = sqrt(255/56) = {float(cbmax):.4f}   [Derived, arbre, A1-A3]")
print(f"  alpha_b^max = c_b^max sqrt(2 w1) = {float(cbmax)*np.sqrt(2*w1):.3f}")
print(f"  TABLE DU CORPUS (Article II) : c_b = 0,30 / 0,60 / 3,02 / 155")
for cb in (0.30,0.60,3.02,155):
    print(f"    c_b = {cb:6.2f} : {'ADMIS' if cb<=float(cbmax) else 'EXCLU par le plongement (arbre)'}")
print("  => le plongement Type I' SELECTIONNE la branche basse du corpus [0,30 ; 0,61]")
print("     et EXCLUT les scenarios 3,02 et 155 — verdict independant de l'EFT.")
print()
print("=== 3. LA COMPOSITION DE Phi (imposee par la fenetre) ===")
w2=np.array([float(c_phi),float(c_chi)]); wn=w2/np.linalg.norm(w2)
wperp=np.array([-wn[1],wn[0]])
print(f"  direction DECOUPLEE de la tension : w_perp = {wperp[0]:+.3f} phi + {wperp[1]:+.3f} chi")
print(f"  etat pur decouple : {wperp[1]**2:.1%} module de volume, {wperp[0]**2:.1%} dilaton")
for s in (0.141,0.286):
    for sgn in (1,-1):
        n=np.cos(np.arcsin(s))*wperp+sgn*s*wn
        pass
vols=[]
for s in (0.141,0.286):
    for sgn in (1,-1):
        n=np.sqrt(1-s*s)*wperp+sgn*s*wn
        vols.append(n[1]**2)
print(f"  avec l'admixture requise (|cos| 0,141-0,286 le long de w) :")
print(f"  FRACTION VOLUME de Phi : {min(vols):.0%} a {max(vols):.0%} ; dilaton : {1-max(vols):.0%} a {1-min(vols):.0%}")
print("  => Phi est PREDIT majoritairement 'module de volume', minoritairement dilaton.")
print()
print("=== 4. LA DISTRIBUTION DE c_b SOUS LA MESURE DE CORDES (le max sans stabilisation) ===")
from scipy.linalg import null_space
rng=np.random.default_rng(20260805); cb=[]
S=np.array([1,-1,2,-2,np.sqrt(2),-np.sqrt(2),1/np.sqrt(2),-1/np.sqrt(2),np.sqrt(7)/2,-np.sqrt(7)/2,2/np.sqrt(7),-2/np.sqrt(7)])
for _ in range(150000):
    A=rng.uniform(0.5,2.0,3)*rng.choice([1,-1],3); E=rng.choice(S,(3,2))
    ns=null_space(E.T)
    if ns.shape[1]!=1: continue
    for sgn in (1,-1):
        u=sgn*ns[:,0]; r=u/A
        if np.any(r<=0): continue
        H=(E.T*u)@E; ev,evec=np.linalg.eigh(H)
        if ev[0]<=0 or ev[0]/ev[1]>0.1: continue
        cb.append(float(cbmax)*abs(evec[:,0]@wn)); break
cb=np.array(cb)
print(f"  vides valides : {len(cb)} ; c_b : mediane {np.median(cb):.3f} ; quartiles [{np.percentile(cb,25):.3f}, {np.percentile(cb,75):.3f}]")
print(f"  P(c_b dans [0,30 ; 0,61]) = {np.mean((cb>=0.30)&(cb<=0.61)):.1%}")
print(f"  P(c_b < 0,30) = {np.mean(cb<0.30):.1%} ; P(c_b > 0,61) = {np.mean(cb>0.61):.1%} ; P(c_b > 2) = {np.mean(cb>2):.1%}")
