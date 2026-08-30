"""Article I — verification de chaque nombre du texte. Boufourou, serie finale."""
import numpy as np
hbar_c=1.97327e-16  # GeV.m
R=8.2e-6; piR=np.pi*R
MPl=2.435e18  # reduite (valeur precise)
M5=(MPl**2/(np.pi*R/hbar_c))**(1/3)  # volume physique de l'orbifold = piR
print(f"pi R = {piR*1e6:.1f} um (texte : 25,8)")
print(f"M5 = (MPl^2/piR)^(1/3) = {M5:.1e} GeV (texte : 3,6e8)")
print(f"rho_Lambda^(1/4) = 2,25 meV (calibration) ; Lambda = C/R^4")
MPl_full=1.2209e19; mKK=hbar_c/R
print(f"echelle des especes (convention litterature, M_Pl NON reduite) :")
print(f"  Lambda_QG = (MPl_full^2 m_KK)^(1/3) = {(MPl_full**2*mKK)**(1/3):.1e} GeV (texte : 1,5e9 ; bande 1e9-1e10 ✓)")
print(f"m_r.R/hbar c = 0,37 => m_r central = {0.37*hbar_c/R*1e12:.1f} meV ; fenetre 5-40 meV")
print(f"portees radion : {hbar_c/40e-12*1e6:.1f}-{hbar_c/5e-12*1e6:.1f} um")
print(f"|1+w| ~ rho_m/(m_r^2 MPl^2) : m_r=9e-12 GeV -> {(1.85e-12)**4/((9e-12)**2*MPl**2):.0e} (~1e-61 attendu, ordre)")
print(f"S ~ MPl^2/m_r = {MPl**2/9e-12:.0e} (>=1e45 ✓)")
print(f"un p.p.m. de deformation : rho ~ 1e-6 m_r^2 MPl^2 = {1e-6*(9e-3)**2*(2.4e27)**2*2.3e-16:.0e} kg/m3 = {1e-6*(9e-3)**2*(2.4e27)**2*2.3e-16/2e17:.0e} x nucleaire")
print(f"vertex mince : (L_g/ell)^2 = {(2e-19/5.8e-6)**2:.0e} (texte : 1e-26)")
print(f"one-scale : T^(1/4)/L_g^-1 = (1/(lambda m))^(1/4) = 1 pour lambda=1/m ✓ (identite)")
print(f"V_stab^(1/4) = sqrt(m_r MPl) : {np.sqrt(5e-12*MPl)/1e3:.1f}-{np.sqrt(40e-12*MPl)/1e3:.1f} TeV (texte 4-11)")
print(f"N_eff : 4 dof, DN=2*(4/7)*x^4<0.3 -> T_d/T_g < {(0.3/(4*4/7))**0.25:.2f} (texte 0,60)")
print(f"        CMB-S4 0.06 -> < {(0.06/(4*4/7))**0.25:.2f} (texte 0,40)")
