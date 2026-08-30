"""Article I §5.3 — le troisieme verrou : la contrainte relique du radion (probleme des modules).
Verifie chaque nombre du paragraphe. Boufourou, serie finale v2.1 (04/08/2026)."""
import numpy as np
MPl=1.22e19            # GeV (masse de Planck)
MPl_red=2.4e18         # GeV (reduite)
gstar=100.0
s0=2891.0              # cm^-3 (entropie aujourd'hui)
rhoc_h2=1.054e-5       # GeV cm^-3
def T_osc(m):          # H = m en radiation
    return np.sqrt(m*MPl/(1.66*np.sqrt(gstar)))
def Omega_h2(m,phi_i):
    rho=0.5*m**2*phi_i**2
    s=(2*np.pi**2/45)*gstar*T_osc(m)**3
    return (rho/s)*s0/rhoc_h2
print("=== Le troisieme verrou : le radion n'a jamais roule de loin ===")
for meV,lab in ((5,"bas de fenetre"),(24,"central"),(40,"haut de fenetre")):
    m=meV*1e-12
    To=T_osc(m); phimax=MPl*np.sqrt(0.12/Omega_h2(m,MPl))
    print(f" m_r = {meV:2d} meV ({lab:14s}) : T_osc = {To/1e3:5.1f} TeV ; "
          f"phi_max = {phimax/MPl:.1e} M_Pl ; dR/R < {phimax/(np.sqrt(6)*MPl):.1e}")
m=24e-12
print(f"\n deplacement planckien : Omega h^2 = {Omega_h2(m,MPl):.1e} (cible 0,12) "
      f"=> surfermeture de {np.log10(Omega_h2(m,MPl)/0.12):.0f} ordres")
print(f" BORNE (texte) : dR/R < {MPl*np.sqrt(0.12/Omega_h2(m,MPl))/(np.sqrt(6)*MPl):.0e} ~ 1e-7")
print("\n=== Remarque structurelle (renvoi Article IV) ===")
MS=7.6e3
print(f" T_osc = sqrt(m M_Pl/1.66 sqrt(g*)) et m = M_S^2/M_Pl (ombre)")
print(f" => T_osc = M_S x sqrt(M_Pl/M_Pl_red)/(1.66 sqrt(g*))^(1/2) ~ {MS*np.sqrt(MPl/MPl_red)/np.sqrt(1.66*np.sqrt(gstar))/1e3:.1f} TeV")
print(f" calcul direct a m = 24 meV : T_osc = {T_osc(24e-12)/1e3:.1f} TeV")
print(f" [accord a un facteur {T_osc(24e-12)/(MS*np.sqrt(MPl/MPl_red)/np.sqrt(1.66*np.sqrt(gstar))):.1f} pres :")
print("  conventions de masse de Planck et g* — c'est un enonce d'ORDRE, declare tel quel]")
print(" => LE RADION COMMENCE A OSCILLER A L'ECHELLE DU MUR : l'ombre revient a sa source.")
