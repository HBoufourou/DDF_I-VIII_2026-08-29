"""EXP-B — G6 : l'amplitude de la Yukawa du capaciteur sigma. Boufourou, 03/08/2026.
Conventions declarees : 5D, [sigma]=GeV^(3/2), [J0]=GeV^(5/2) (BC sigma'=J0/2),
[g]=GeV^(1/2) (couplage g.Phi^2.sigma). k_corpus = 8.16e-22 GeV^2 (pente E/longueur)."""
import numpy as np
GeV_inv_m=1.97327e-16   # m par GeV^-1
piR=25.8e-6/GeV_inv_m   # GeV^-1
rhoL=(2.25e-12)**4      # GeV^4 (2.25 meV)
M5=3.5e8                # GeV
k=8.161e-22             # GeV^2

print("=== ETAPE 1 : la matiere source-t-elle sigma ? (le coeur de G6) ===")
print("  Dans le modele minimal {sigma, tadpoles J0/J1, g.Phi^2.sigma}, le tadpole")
print("  J0.sigma.delta(y) couple sigma a l'OPERATEUR UNITE de la brane (un terme de")
print("  tension), PAS a la densite de matiere. Un atome ne porte aucune 'charge sigma'.")
print("  ==> a l'ordre des arbres, une masse d'epreuve NE SOURCE PAS sigma.")
print("      FORCE ENTRE DEUX MASSES DE LABO : ZERO. AUCUNE raie sigma dans (P1b).\n")
print("=== ETAPE 2 : les canaux induits (a chiffrer avant de dire 'zero') ===")
# (a) melange sigma-radion : le fond de sigma depend de R => melange de masse
V_sig=(30e-3*1e-9)**4   # (30 meV)^4 en GeV^4
V_stab=(4.6e3)**4       # (4.6 TeV)^4 en GeV^4
theta=V_sig/V_stab
print(f"  (a) melange sigma-radion : theta ~ V_sigma/V_stab = {theta:.0e}")
print(f"      le radion couple a 1/M_Pl ; sigma herite ~ theta/M_Pl =>")
print(f"      rapport de force a la gravite ~ theta^2 = {theta**2:.0e}  ->  RIEN (1e-113).")
# (b) fond de Phi dans un halo : sigma source par g<Phi^2> => force noir-noir
print(f"  (b) g.<Phi^2> dans un halo source sigma : force ENTRE CONCENTRATIONS SOMBRES")
print(f"      uniquement — invisible au labo, correction auto-coherente du puits ailleurs.")
print(f"  (c) couplage a un operateur SM (sigma.|H|^2, sigma.T) : ABSENT du modele")
print(f"      minimal ; M4 exige une justification ecrite avant tout ajout. DECLARE.\n")
print("=== ETAPE 3 : les VRAIES contraintes residuelles du secteur sigma ===")
# energie du fond : champ (1/2)(sigma')^2 sur l'intervalle + energie des sources
u_field=0.5*(0.5)**2    # (1/2)(J0/2)^2 par unite J0^2
E_field=u_field*piR     # x J0^2  -> GeV^-1 * ... => rho_4D = (J0^2/8) piR  [GeV^4 si J0^2 GeV^5 x GeV^-1]
# rho_sigma = (piR/8) J0^2  <= rho_Lambda  => J0max
J0max=np.sqrt(8*rhoL/piR)
print(f"  (i) budget Lambda : rho_sigma = (piR/8).J0^2 <= rho_Lambda = {rhoL:.1e} GeV^4")
print(f"      => J0 <= {J0max:.1e} GeV^(5/2)")
g_min=2*k/J0max
print(f"      et k = g.J0/2 fixe => g >= 2k/J0max = {g_min:.1e} GeV^(1/2)")
print(f"      naturalite : sqrt(M5) = {np.sqrt(M5):.1e} GeV^(1/2)  =>  g/sqrt(M5) >= {g_min/np.sqrt(M5):.0f}")
print(f"      ==> TENSION DE NATURALITE ~ O(10^3) A DECLARER : le capaciteur demande")
print(f"          soit un couplage ~2000x super-naturel, soit une annulation partielle")
print(f"          de son energie de fond (territoire S-3, herite, pas cree ici).")
# (ii) terme lineaire en R du potentiel inter-branes
E_tot=-3/8  # x J0^2 piR (champ + sources, capaciteur a charges opposees)
print(f"  (ii) le secteur sigma ajoute au potentiel inter-branes un terme LINEAIRE en R")
print(f"       (force constante entre les branes), d'amplitude <= rho_Lambda ~ (2 meV)^4 :")
print(f"       {V_stab/rhoL:.0e} fois sous le stabilisateur TeV — il ne stabilise ni ne")
print(f"       destabilise rien de mesurable (coherent avec A2/A3).\n")
print("=== VERDICT G6 ===")
print("  L'amplitude de la Yukawa du sigma entre masses d'epreuve est NULLE a l'ordre")
print("  des arbres et < 1e-100 x gravite par les canaux induits : LA RAIE SIGMA EST")
print("  INVISIBLE — le spectre (P1b) est INCHANGE. Le brouillon peut l'ecrire ainsi :")
print("  'the source field of the well leaves no fifth-force trace; its only residues")
print("   are a Lambda-budget bound on the tadpole (J0 <= 4e-29 GeV^(5/2), implying a")
print("   coupling ~10^3 above sqrt(M5), declared) and an R-linear inter-brane term")
print("   bounded by the dark-energy density.'  ==> LE BLOQUANT DU BROUILLON EST LEVE.")
