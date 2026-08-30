"""LVS DE REFERENCE — implementation validee, minimum VERIFIE. 07/08/2026.
Corrige deux erreurs : (1) l'alpha' correction porte xi_hat = xi/g_s^{3/2}, pas xi ;
(2) ma relation tau_s ~ xi^{2/3}/g_s etait approximative (facteur 2^{2/3})."""
import numpy as np
from scipy.optimize import brentq, minimize
MPl=2.435e18; meV=1e-12

def V_LVS(lnV, tau_s, xi, g_s, a_s, A_s, W0, lam=1.0):
    """Potentiel LVS standard. xi_hat = xi/g_s^{3/2} (correction alpha' habillee du dilaton)."""
    V=np.exp(lnV); xih=xi/g_s**1.5
    t1=(8/3)*(a_s*A_s)**2*np.sqrt(tau_s)*np.exp(-2*a_s*tau_s)/(lam*V)
    t2=-4*a_s*A_s*W0*tau_s*np.exp(-a_s*tau_s)/V**2
    t3=3*xih*W0**2/(4*V**3)
    return t1+t2+t3

def minimum_LVS(xi,g_s,N,A_s=1.0,W0=1e-9,lam=1.0):
    a_s=2*np.pi/N
    # depart : relations analytiques CORRECTES
    xih=xi/g_s**1.5
    tau0=(xih/(2*lam))**(2/3)
    lnV0=np.log(3*np.sqrt(tau0)*lam*W0*np.exp(a_s*tau0)/(4*a_s*A_s))
    f=lambda z: V_LVS(z[0],z[1],xi,g_s,a_s,A_s,W0,lam)
    r=minimize(f,[lnV0,tau0],method='Nelder-Mead',
               options={'xatol':1e-12,'fatol':1e-40,'maxiter':20000,'maxfev':20000})
    lnV,ts=r.x
    # audit sans dimension
    h=1e-5; g=np.array([(f([lnV+h,ts])-f([lnV-h,ts]))/(2*h),(f([lnV,ts+h*ts])-f([lnV,ts-h*ts]))/(2*h*ts)])
    sc=np.array([1.0,ts]); Vm=abs(r.fun) if r.fun!=0 else 1
    H=np.zeros((2,2)); hh=[1e-4,1e-4*ts]
    for i in range(2):
        for j in range(2):
            e1=np.zeros(2); e1[i]=hh[i]; e2=np.zeros(2); e2[j]=hh[j]
            H[i,j]=(f(r.x+e1+e2)-f(r.x+e1-e2)-f(r.x-e1+e2)+f(r.x-e1-e2))/(4*hh[i]*hh[j])
    Hl=H*np.outer(sc,sc)/Vm; ev=np.linalg.eigvalsh((Hl+Hl.T)/2)
    return dict(V=np.exp(lnV),tau_s=ts,a_s_tau_s=a_s*ts,Vmin=r.fun,
                grad_log=np.linalg.norm(g*sc/Vm),eigs=ev,
                minimum=bool(np.linalg.norm(g*sc/Vm)<1e-3 and ev.min()>0),
                m32_meV=W0*MPl/np.exp(lnV)/meV)

print("="*76); print("A. LES DEUX CORRECTIONS QUI EXPLIQUENT TES TROIS RESULTATS CONTRADICTOIRES")
print("="*76)
print("  (1) la correction alpha' porte xi_hat = xi/g_s^{3/2}, PAS xi.")
for xi,gs in ((0.5,0.02),(1.0,0.02),(2.0,0.05)):
    print(f"      xi={xi:4.1f}, g_s={gs:5.3f}  ->  xi_hat = {xi/gs**1.5:8.1f}   [ton xi=10990 d'hier etait un xi_HAT]")
print("  (2) la relation exacte est tau_s^{3/2} = xi_hat/2, soit tau_s = xi^{2/3}/(2^{2/3} g_s).")
print("      Ma formule tau_s ~ xi^{2/3}/g_s etait trop grande d'un facteur 2^{2/3} = 1,587.")
for xi,gs in ((0.5,0.02),):
    print(f"      xi=0,5 g_s=0,02 : ma formule -> {xi**(2/3)/gs:.1f} ; EXACTE -> {(xi/(2*gs**1.5))**(2/3):.1f}"
          f"   (ton sympy donnait 18,2 : c'est LUI qui avait raison)")

print("\n"+"="*76); print("B. LE POINT DE PARAMETRES CORRECT (resolu, pas ajuste)")
print("="*76)
def xi_pour(target_atau,g_s,N,lam=1.0):
    a_s=2*np.pi/N; tau=target_atau/a_s
    return 2*lam*g_s**1.5*tau**1.5
print(f"  cible : a_s tau_s = 66,3  (pour m_3/2 ~ 24 meV)")
print(f"  {'N':>3} {'g_s':>7} {'tau_s':>8} {'xi requis':>10} {'|chi| requis':>13} {'realiste ?':>11}")
for N in (2,3,4,6):
    for gs in (0.01,0.02,0.05,0.1):
        a_s=2*np.pi/N; tau=66.3/a_s; xi=xi_pour(66.3,gs,N)
        chi=2*xi*(2*np.pi)**3/1.2021
        ok = "OUI" if (0.05<=xi<=2.4 and chi<=960) else "non"
        if ok=="OUI": print(f"  {N:3d} {gs:7.3f} {tau:8.2f} {xi:10.3f} {chi:13.0f} {ok:>11}")
print("\n"+"="*76); print("C. VERIFICATION DIRECTE : minimum trouve et AUDITE")
print("="*76)
for N,gs in ((3,0.02),(2,0.01),(4,0.05)):
    xi=xi_pour(66.3,gs,N)
    if not(0.05<=xi<=2.4): continue
    r=minimum_LVS(xi,gs,N,W0=1.44e-9)
    print(f"  N={N}, g_s={gs:5.3f}, xi={xi:.3f}  ->  V = {r['V']:.3e} ; tau_s = {r['tau_s']:.2f} ;"
          f" a_s tau_s = {r['a_s_tau_s']:.1f}")
    print(f"       m_3/2 = {r['m32_meV']:.1f} meV ; |grad log| = {r['grad_log']:.1e} ;"
          f" val. propres = {np.array2string(r['eigs'],precision=2)} -> MINIMUM : {r['minimum']}")
