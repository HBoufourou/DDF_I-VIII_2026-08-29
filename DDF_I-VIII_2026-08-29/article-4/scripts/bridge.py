"""Reproduces every number of the bridge paper."""
import numpy as np
# p from the response chain
print("response chain: |grad|~sqrt(alpha) then recoil x alpha  ->  a_D ~ alpha^(3/2)  =>  p=3/2")
deficit=6e-11/2e-22
print(f"deficit on g_D = {deficit:.1e} ; coupling target = deficit^(2/3) = {deficit**(2/3):.1e}")
# condensate occupation (obstruction I closure, from companion)
m=0.025; n1=1.7e-6; v=1e-3
N0=n1*(2*np.pi/(m*v))**3
print(f"ground-level occupation n1*lambda_dB^3 = {N0:.1e}  (threshold 2.61: {np.log10(N0/2.61):.0f} orders above)")
# discriminant table
b1=1.019
for lab,u in [("posited-field",0),("gravitational u=1",1.0),("gravitational u=2.3",2.3),("gravitational u=4.44",4.44)]:
    w=1/b1 if u==0 else np.exp(-(4/3)*u**1.5)
    print(f"{lab:22s}: |psi1(0)|^2={w:.2e}  g5 required={1/w:.1e}")
umax=(0.75*np.log(100))**(2/3)
print(f"structural bound (g5<=100): u <= {umax:.2f}  ->  <= {(2/(3*np.pi))*umax**1.5:.2f} even level per branch")

# --- juge astrophysique (SN1987A) et borne de portail ---
T_core=30e6; MPl_full=1.22e28; yN=0.3*0.939/246
print(f"SN1987A : suppression (T/M_Pl)^2 = {(T_core/MPl_full)**2:.1e} -> passe trivialement (couplage gravitationnel)")
print(f"  portail : SN1987A donnerait theta < {1e-10/yN:.1e} ; les forces courtes donnent theta < 4e-16 (8 ordres plus fort)")

# --- verification symbolique du paragraphe covariant (02/08 soir) ---
import sympy as sp
cb,w1,R,M5=sp.symbols('c_b w_1 R M_5',positive=True)
g1=cb*sp.sqrt(w1/(sp.pi*R))/M5**sp.Rational(3,2)
alpha=sp.simplify(g1*sp.sqrt(M5**3*2*sp.pi*R))
assert sp.simplify(alpha-cb*sp.sqrt(2*w1))==0
print("operateur covariant : alpha = c_b sqrt(2 w1) — identite verifiee")
