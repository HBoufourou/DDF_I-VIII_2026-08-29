"""V2 (faille A) + V3 (faille B) — resolution structurelle + la contrainte quantitative
qui en sort. 05/08/2026. Chaque nombre de la note sort d'ici."""
import numpy as np
MPl=2.435e18; eps=24.0  # meV, unite d'Airy = premier niveau (corpus)
r=[1.0,2.294,3.190,4.031]  # rapports de la tour (corpus, valides)
print("=== V3 : LES DEUX MECANISMES DE MASSE — la combinaison, pas le choix ===")
print("  spectre 4D d'un champ 5D avec masse molle m0 (uniforme) + puits :")
print("  m_n^2 = m0^2 + E_n^2   =>   les RAPPORTS de la tour se deforment si m0 ~ eps\n")
print(f"  {'m0/eps':>7} {'m2/m1':>7} {'m3/m1':>7} {'m4/m1':>7}   (corpus : 2,294  3,190  4,031)")
for x in (0.0,0.1,0.22,0.5,1.0):
    m=[np.sqrt(x**2+ri**2) for ri in r]
    print(f"  {x:7.2f} {m[1]/m[0]:7.3f} {m[2]/m[0]:7.3f} {m[3]/m[0]:7.3f}")
print("\n  tolerance sur r2 -> borne sur m0 (developpement + exact) :")
for tol in (0.01,0.02,0.05):
    # exact : resoudre r2(m0)=2.294(1-tol)
    from scipy.optimize import brentq
    f=lambda x: np.sqrt(x*x+r[1]**2)/np.sqrt(x*x+1)-r[1]*(1-tol)
    x0=brentq(f,0,3)
    print(f"    |dr2/r2| < {tol:4.0%}  =>  m0 < {x0:.2f} eps = {x0*eps:.1f} meV")
print("\n=== LA PROTECTION NO-SCALE, ET LA CONTRAINTE CROISEE QUI EN SORT ===")
print("  K = -3 ln(T+Tb) est NO-SCALE : la masse molle du module de VOLUME s'annule a")
print("  l'arbre [Cremmer-Ferrara-Kounnas-Nanopoulos 1983] ; le DILATON n'est pas protege")
print("  => m0(Phi) ~ |sin theta_d| x m_3/2   (theta_d = fraction dilaton de Phi)")
print("  => la tour n'est integre que si sin^2(theta_d) < (m0_max/m_3/2)^2\n")
print(f"  {'M_S (TeV)':>10} {'m_3/2 (meV)':>12} {'fraction dilaton max (tol 2%)':>30}")
from scipy.optimize import brentq
f=lambda x: np.sqrt(x*x+r[1]**2)/np.sqrt(x*x+1)-r[1]*0.98
m0max=brentq(f,0,3)*eps
for MS in (4,7.6,10):
    m32=(MS*1e3)**2/MPl*1e12
    fd=min(1,(m0max/m32)**2)
    print(f"  {MS:10.1f} {m32:12.1f} {fd:29.0%}")
print(f"\n  => a M_S = 4 TeV : TOUTE la plage de composition (dilaton 9-60 %) survit ;")
print(f"     a M_S = 7,6-10 TeV : la fraction dilaton est ECRASEE a quelques % ;")
print(f"     la composition d'E3-FINAL (volume 40-91 %) et l'integrite de la tour")
print(f"     se COMPATIBILISENT en poussant Phi vers le volume — la meme direction")
print(f"     que la protection no-scale. TROIS logiques, une conclusion : Phi ~ volume.")
