"""Article IV — verification de chaque nombre du texte (C5b + C7)."""
import numpy as np
MPl=2.4e18
print("== la relation d'ombre ==")
for MS in (3.5e3,4.6e3,7.6e3,9.8e3):
    print(f"  M_S = {MS/1e3:4.1f} TeV -> mu = M_S^2/MPl = {MS**2/MPl*1e12:5.1f} meV")
print(f"  1/R = 24 meV -> M_S = sqrt(24e-12*MPl) = {np.sqrt(24e-12*MPl)/1e3:.1f} TeV (texte 7,6)")
print(f"  fenetre radion 5-40 meV -> M_S = {np.sqrt(5e-12*MPl)/1e3:.1f}-{np.sqrt(40e-12*MPl)/1e3:.1f} TeV")
print("== la convergence ==")
T14=9.8  # convention piR, M5 = 3,568e8 (u_derive.py)
print(f"  T^(1/4) (Article II) = {T14} TeV vs knot lifte 7,6 TeV : accord a {abs(T14-7.6)/7.6*100:.0f}% (texte ~30%)")
print(f"  [la precision qu'autorisent les coefficients O(1) declares : le 1/6 de la pente, le O(1) de l'ombre]")
print("== le spectre de l'ombre ==")
mu=np.sqrt(24e-12*MPl)**2/MPl*1e12  # = 24.0... recompute proprement :
mu=(7.6e3)**2/MPl*1e12
print(f"  mu (ancre a M_S=7,6) = {mu:.1f} meV (texte 24,1)")
for nom,m in (("m_Phi",24),("m_r bas",5),("m_r haut",40),("m_sigma max",4.3),("m_3/2",mu/np.sqrt(3))):
    print(f"  {nom:12s} = {m:5.1f} meV -> c_i = {m/mu:.2f}")
print("== le cross-check de boucle ==")
print(f"  Lambda^(1/4)/mu = 2.25/{mu:.1f} = {2.25/mu:.3f} vs 1/4pi = {1/(4*np.pi):.3f}")
print("== le triangle ==")
print(f"  balance : m_r mesure -> M_S = sqrt(m_r MPl) a la largeur de c_r pres")
print(f"  beta/H in [10,100] conditionnel au 1er ordre (S3/T ~ 140 a f ~ 6 TeV)")
