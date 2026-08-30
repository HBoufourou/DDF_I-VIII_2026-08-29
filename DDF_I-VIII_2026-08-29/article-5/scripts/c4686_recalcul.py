import numpy as np, sympy as sp
np.set_printoptions(precision=5, suppress=False, linewidth=150)
kap=np.zeros((5,5,5))
for (i,j,k),v in {(0,3,4):2,(1,1,1):2,(1,1,3):-2,(1,3,3):2,(2,2,2):2,(2,2,4):-2,
                  (2,4,4):2,(3,3,3):-2,(3,3,4):-2,(4,4,4):-2}.items():
    for p in set([(i,j,k),(i,k,j),(j,i,k),(j,k,i),(k,i,j),(k,j,i)]): kap[p]=v
def V_of(t):  return np.einsum('ijk,i,j,k',kap,t,t,t)/6.
def tau_of(t):return np.einsum('ijk,j,k',kap,t,t)/2.
def M_of(t):  return np.einsum('ijk,k',kap,t)              # dtau_i/dt_j
print("="*78); print("1. LE POINT LVS — resolu proprement (tau_s > 0 strict)")
print("="*78)
from scipy.optimize import fsolve
Vt=1e6; rs=8.6411e-12; tw=0.5     # tau_s/V^(2/3) et tau_W/V^(2/3)
def eqs(t):
    tau=tau_of(t)
    return [V_of(t)-Vt, tau[1]-rs*Vt**(2/3), tau[3]-tw*Vt**(2/3),
            t[2]-0.40*Vt**(1/3), t[4]-0.62*Vt**(1/3)]
t0=np.array([135.8,94.15,40.,94.15,62.])
sol,info,ier,msg=fsolve(eqs,t0,full_output=True,xtol=1e-14)
tau=tau_of(sol); V=V_of(sol)
print(f"  t   = {sol}")
print(f"  tau = {tau}")
print(f"  V = {V:.6e} ; tau_s = {tau[1]:.4e} ; tau_s/V^(2/3) = {tau[1]/V**(2/3):.3e}  (cible 8,64e-12)")
print(f"  convergence : {ier==1} ; tau_s STRICTEMENT positif : {tau[1]>0}")
print("\n"+"="*78); print("2. LA METRIQUE DE KAHLER DANS LA BASE DES MODULES tau (pas t)")
print("="*78)
M=M_of(sol); Minv=np.linalg.inv(M)
# dK/dtau_i = -t_i/V  ;  d2K/dtau_i dtau_j = -(dt_i/dtau_j)/V + t_i t_j /V^2
G=(-Minv/V + np.outer(sol,sol)/V**2)/1.0
G=0.5*(G+G.T)
ev,U=np.linalg.eigh(G)
print(f"  eig(G_tau) = {ev}")
print(f"  definie positive : {np.all(ev>0)} ; conditionnement = {ev.max()/ev.min():.2e}")
print("  (dans SA base t, la petite valeur propre etait un artefact de coordonnees :")
print("   en base tau le blow-up qui retrecit donne une GRANDE valeur propre, pas une petite)")
print("\n"+"="*78); print("3. DIRECTIONS FIXEES, DIRECTIONS PLATES")
print("="*78)
gV=sol/(2*V)                      # d(ln V)/d tau_i = t_i/(2V)
eS=np.eye(5)[1]                   # tau_s (dP7 LVS)  fixe par le LVS
eD=np.eye(5)[3]                   # Wilson D-terme   fixe par le D-terme
fixes=[gV,eS,eD]
B=[]
for e in np.eye(5):
    v=e.copy()
    for f in fixes+B: v=v-(v@G@f)/(f@G@f)*f
    if np.sqrt(abs(v@G@v))>1e-10*np.sqrt(abs(e@G@e)): B.append(v/np.sqrt(v@G@v))
print(f"  directions plates trouvees : {len(B)}   (attendu 5-2-1 = 2)")
for i,b in enumerate(B): print(f"    n{i+1} = {b}")
print("\n"+"="*78); print("4. LES PROJECTIONS p1, p2 — enroulement sur le Wilson (index 4)")
print("="*78)
wg=1.6956
for D,lab in ((4,"Wilson enroulement (prime 8)"),(0,"K3 (prime 1)"),(2,"dP7 plate (prime 6)")):
    cov=np.eye(5)[D]/tau[D]                       # d(ln tau_D)/d tau
    nrm=np.sqrt(cov@np.linalg.inv(G)@cov)
    ps=[wg*(cov@b)/nrm for b in B]
    R=np.sqrt(sum(p*p for p in ps))
    tt=np.linspace(0,2*np.pi,400000); cb=np.abs(np.cos(tt)*ps[0]+np.sin(tt)*ps[1]) if len(ps)>1 else np.abs(ps[0])*np.ones(1)
    fr=float(np.mean((cb>=0.30)&(cb<=0.61)))
    print(f"  {lab:30s} p = {np.round(ps,4)}  |p| = {R:.4f}  fenetre atteignable : {R>0.30}  fraction {fr:.1%}")
