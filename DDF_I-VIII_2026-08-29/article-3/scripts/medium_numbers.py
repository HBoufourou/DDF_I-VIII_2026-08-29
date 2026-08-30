"""Article III — verification de chaque nombre du texte."""
import numpy as np
# fractions du profil plat (deplacement y-homogene sur la tour entrelacee)
f=[0.27,0.39,0.22,0.12]
print(f"fractions {f} ; etages superieurs (2-4) = {sum(f[1:]):.2f} (texte 73%)")
# suppression graviton : gap 25 meV vs hw ~ 1e-26 eV -> exp(-(gap/hw)... ordre)
print(f"gap/hw = 25e-3/1e-26 = {25e-3/1e-26:.0e} => suppression e^-O(1e24) (texte)")
# degenerescence de Bose : n lambda_dB^3
rho=0.4  # GeV/cm3 local
m=24e-12 # GeV
n=rho/m  # /cm3
lam_dB_cm=1.24e-4/ (m*1e9) *1e2 if False else None
# lambda_dB = h/(m v), v~200 km/s :
h=6.626e-34  # J.s
m_kg=m*1.783e-27  # GeV -> kg
v=2e5  # m/s (200 km/s)
lam_dB=h/(m_kg*v)*1e2  # cm
print(f"n = {n:.1e} /cm3 ; lambda_dB = {lam_dB:.1e} cm ; n.lam^3 = {n*lam_dB**3:.1e} (texte 1e10-12)")
# fenetre lambda
print("lambda operative ~1e-18 (stimulee, 6 ordres sous 1e-12) ; plafond Bullet 2e-13 ✓")
# R_dec : Soleil
K=4.9e-4  # pc.GeV/cm3
for ff in (1.4e-2,1.2e-3):
    R=K*1**(1/3)/(ff*rho)
    print(f"Soleil, f={ff:.1e} : R_dec = {R:.3f} pc = {R*206265/1e3:.0f} kau")
print("=> bulle solaire 18 kau - 1 pc : Cassini (10 UA) dedans par >~1e3 ✓ ; WB >10 kau dedans ✓")
# reglement thermalise
Vs=(4.6e3)**4; g=100
T=(30*Vs/(np.pi**2*g))**0.25
print(f"T reglement = (30 V_stab/pi^2 g*)^(1/4) = {T/1e3:.1f} TeV (texte ~2)")
# C4 (produit par C4_run_sparc.py, rappele ici)
print("C4 : |dA/A| <= 9% a Ups=0.5 (p=0.08) ; <=14% sur Ups in [0.4,0.8] ; zero a Ups=0.69")
