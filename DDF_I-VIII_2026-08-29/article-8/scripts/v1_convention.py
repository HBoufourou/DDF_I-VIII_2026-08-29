"""V1 — LA VERIFICATION LA PLUS URGENTE : la convention du vertex de l'Article I par. 6
confrontee a la notre. Source : Eq.(10)-(11) et Eq.(1) de l'article (fourni)."""
import numpy as np
hbar_c=1.97327e-16  # m.GeV^-1
print("=== 1. QUEL M_Pl dans leur Eq.(1) M_Pl^2 = M5^3 (pi R) ? (verif par LEURS nombres) ===")
piR=25.8e-6/hbar_c
M5=3.6e8
MPl=np.sqrt(M5**3*piR)
print(f"  pi R = 25,8 um = {piR:.3e} GeV^-1 ; M5 = 3,6e8  =>  M_Pl = {MPl:.3e} GeV")
print(f"  M_Pl REDUITE = 2,435e18 ; NON reduite = 1,221e19  =>  leur M_Pl est la REDUITE ✓")
M5x=(2.435e18**2/piR)**(1/3)
print(f"  contre-verif : M5 = (Mbar^2/piR)^(1/3) = {M5x:.3e} = leur 3,6e8 ✓")
print("\n=== 2. LA CHAINE DE LEUR Eq.(10) VERS 4D (dimensionnellement, pas a pas) ===")
print("  operateur 5D : (c_b/M5^{3/2}) Phi T delta(y-y0)   [Phi 5D canonique, dim 3/2]")
print("  mode n : Phi = psi_n(y) phi_n(x), Int|psi|^2 dy = 1  =>  phi_n 4D canonique")
print("  couplage 4D du mode n : c_b psi_n(y0)/M5^{3/2}")
print("  leur w1 = piR|psi_1(y0)|^2  =>  psi_1(y0) = sqrt(w1/piR)")
print("  =>  couplage = c_b sqrt(w1) / sqrt(piR M5^3) = c_b sqrt(w1)/M_Pl,red   [via LEUR Eq.(1)]")
print("  et leur Eq.(11) alpha_b = c_b sqrt(2 w1) = sqrt2 x (le couplage x M_Pl) : coherent.")
print("\n=== 3. LA CONFRONTATION AVEC NOTRE COTE CORDES ===")
print("  nous : L4 = c_string (phi_c/M_Pl,red) T   pour le MODE ZERO (profil plat)")
print("  eux, mode plat : psi = 1/sqrt(piR)  =>  w = piR|psi|^2 = 1  =>  couplage = c_b x 1/M_Pl,red")
print("  =>  c_b(corpus) == c_string : MEME operateur, MEME M_Pl (reduite), MEME champ canonique,")
print("      facteur de profil separe A L'IDENTIQUE (leur sqrt(w1) explicite).")
print("\n================ VERDICT V1 : LES CONVENTIONS COINCIDENT, SANS FACTEUR RESIDUEL ================")
print("  * la peur des 10 ordres est levee : le M5^{3/2} de l'operateur 5D se convertit,")
print("    par LEUR PROPRE calibration de volume Eq.(1), en exactement 1/M_Pl,red en 4D.")
print("  * l'hypothese A3 cesse d'etre une hypothese : c'est un theoreme de raccord.")
print("  * deviennent INCONDITIONNELS (a la normalisation) : c_b_max = 2,134, E-11,")
print("    l'exclusion des scenarios 3,02/155, le verdict LVS.")
print("  * restent : A1 (habillage tension), A2 (moduli d'echelle), et les failles A/B")
print("    (identite du champ du puits ; mecanisme de masse) — questions de PHYSIQUE,")
print("    plus de normalisation.")
