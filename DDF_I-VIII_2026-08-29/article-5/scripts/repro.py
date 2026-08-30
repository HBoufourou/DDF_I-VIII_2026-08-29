"""REPRODUCTION COMPLETE — tous les nombres-cles du dossier CDD-cordes.
Chaque bloc ecrit un JSON dans 06_preuves/. Executable seul : python3 repro.py"""
import numpy as np, sympy as sp, json, os
OUT="/mnt/user-data/outputs/cdd-cordes/06_preuves"; os.makedirs(OUT,exist_ok=True)
def save(name,d): json.dump(d,open(f"{OUT}/{name}.json","w"),indent=1,default=float); print(f"  -> {name}.json")
sq=sp.sqrt

# ============================================================ 1. LES TROIS VECTEURS
print("[1] vecteurs de couplage et identites de norme")
a1=1/(4*sq(7)); b1=-sq(7)/4; a2=sq(sp.Rational(5,28)); b2=-1/sq(35)
assert sp.simplify(7*a1+b1)==0 and sp.simplify(2*a2+5*b2)==0     # frame d'Einstein
w=[sp.simplify(x*sq(2)) for x in (sp.Rational(5,4), 4*a2+5*b2, 9*a1)]           # tension D8
v=[sp.simplify(x*sq(2)) for x in (sp.Rational(5,2), 4*a2+5*b2, 9*a1+b1)]        # flux F0^2
u=[sp.simplify(-2*x*sq(2)) for x in (sp.Rational(-3,2), 5*b2, 5*a1+b1)]         # mur de domaine
nw=sp.simplify(sum(x**2 for x in w)); nv=sp.simplify(sum(x**2 for x in v)); nu=sp.simplify(sum(x**2 for x in u))
wv=sp.simplify(sum(a*b for a,b in zip(w,v))); wu=sp.simplify(sum(a*b for a,b in zip(w,u)))
save("01_vecteurs",dict(w=[float(x) for x in w],v=[float(x) for x in v],u=[float(x) for x in u],
   norme_w2=int(nw),norme_v2=int(nv),norme_u2=int(nu),w_dot_v=int(wv),w_dot_u=int(wu),
   angle_wv_deg=float(sp.deg(sp.acos(wv/sq(nw*nv)))),angle_wu_deg=float(sp.deg(sp.acos(wu/sq(nw*nu)))),
   controle_GM="c_phi = (p-3)/(2 sqrt2), identite en p"))

# ============================================================ 2. BORNE ET FENETRES
print("[2] borne c_b^max, E-11, bandes bimodales")
W=np.array([float(x) for x in w]); cbmax=float(np.hypot(W[0],W[1])); psiw=np.degrees(np.arctan2(W[1],W[0]))
bandes=[]
for lo,hi in ((psiw+73.4,psiw+81.9),(psiw-81.9,psiw-73.4)):
    bandes.append(sorted([np.sin(np.radians(lo))**2, np.sin(np.radians(hi))**2]))
save("02_borne_fenetres",dict(c_b_max=cbmax, c_b_max_exact="sqrt(255/56)", psi_w_deg=float(psiw),
   E11_deg=[73.4,81.9], bandes_p_volume=sorted(bandes), fenetre_c_b=[0.30,0.61],
   benchmarks_exclus=[3.02,155.0], w1=5.5, alpha_b_sur_c_b=float(np.sqrt(2*5.5))))

# ============================================================ 3. VIDE LVS DE REFERENCE
print("[3] vide LVS de reference")
MPl=2.435e18; meV=1e-12
def V_LVS(lnV,ts,xi,gs,a_s,A_s,W0,lam=1.0):
    V=np.exp(lnV); xih=xi/gs**1.5
    return ((8/3)*(a_s*A_s)**2*np.sqrt(ts)*np.exp(-2*a_s*ts)/(lam*V)
            -4*a_s*A_s*W0*ts*np.exp(-a_s*ts)/V**2 + 3*xih*W0**2/(4*V**3))
from scipy.optimize import minimize
vides=[]
for N,gs,W0 in ((2,0.010,1.32e-9),(2,0.020,1.30e-9),(3,0.020,2.60e-9)):
    a_s=2*np.pi/N; tau=66.3/a_s; xi=2*gs**1.5*tau**1.5; xih=xi/gs**1.5; tau0=(xih/2)**(2/3)
    f=lambda z: V_LVS(z[0],z[1],xi,gs,a_s,1.0,W0)/1e-90
    z0=[np.log(3*np.sqrt(tau0)*W0*np.exp(a_s*tau0)/(4*a_s)),tau0]
    r=z0
    for m in ('Nelder-Mead','Powell','Nelder-Mead'):
        r=minimize(f,r,method=m,options={'maxiter':60000,'maxfev':60000}).x
    V=np.exp(r[0]); vides.append(dict(N=N,g_s=gs,xi=float(xi),chi=float(2*xi*(2*np.pi)**3/1.2021),
        volume=float(V),tau_s=float(r[1]),a_s_tau_s=float(a_s*r[1]),m32_meV=float(W0*MPl/V/meV),
        M_S_TeV=float(np.sqrt(W0*MPl/V*MPl)/1e3)))
save("03_vide_LVS",dict(vides=vides,relation="tau_s^{3/2} = xi_hat/2 avec xi_hat = xi/g_s^{3/2}",
   sensibilite="d ln m32 / d ln g_s = + a_s tau_s = 66,7"))

# ============================================================ 4. LES c_b DERIVES
print("[4] les c_b derives (le motif)")
pV=float(7*np.sqrt(6)/12); pF=float(np.sqrt(5/6)); pU=float(11*np.sqrt(6)/12)
motif=[dict(construction="LVS fibre (fibre pure)",c_b=pF,forme="sqrt(5/6)"),
       dict(construction="LVS vanille (volume global)",c_b=pV,forme="7 sqrt6/12"),
       dict(construction="mur de domaine",c_b=pU,forme="11 sqrt6/12"),
       dict(construction="variante w+u",c_b=float(17*np.sqrt(13)/26),forme="17 sqrt13/26"),
       dict(construction="source alignee sur w",c_b=float(np.sqrt(6)),forme="sqrt6")]
for m in motif: m["alpha_b"]=m["c_b"]*float(np.sqrt(2*5.5)); m["dans_fenetre"]=bool(0.30<=m["c_b"]<=0.61)
save("04_motif",dict(motif=motif,fenetre=[0.30,0.61],
   theoreme_negatif="aucun etat propre PUR ne satisfait E-11 ; il faut un melange",
   croisement_conditionnel=dict(tau_f_etoile="sqrt(2 A V)",c_b_antisym=float((pV-pF)/np.sqrt(2)),
       c_b_sym=float((pV+pF)/np.sqrt(2)))))

# ============================================================ 5. TADPOLES / IIA MASSIVE
print("[5] verdict des tadpoles")
hbar_c=1.97327e-16; R=8.2e-6; gs=0.010
tad=[]
for Ms in (0.8e8,2.0e8,3.6e8):
    ls=hbar_c/Ms; r8=gs**2*R/(2*ls)
    tad.append(dict(M_s_GeV=Ms,l_s_m=ls,R_sur_ls=R/ls,r_F0_8=r8,r_F0_1=r8/64))
save("05_tadpoles",dict(formule="r = c^2 g_s^2 R / (2 l_s)",points=tad,cible_r=[0.149,0.419],
   verdict="IIA massive incompatible avec un intervalle mesoscopique (13-16 ordres)",
   configuration="F0 = 0 obligatoire => annulation locale => 8 D8 + O8 a chaque bord",
   consequence="la paroi lointaine porte un secteur de jauge cache, contraint par Delta N_eff"))

# ============================================================ 6. CANDIDAT Id 4686
print("[6] candidat Id 4686 : cone, metrique, directions plates, projections")
kap=np.zeros((5,5,5))
for (i,j,k),val in {(0,3,4):2,(1,1,1):2,(1,1,3):-2,(1,3,3):2,(2,2,2):2,(2,2,4):-2,
                    (2,4,4):2,(3,3,3):-2,(3,3,4):-2,(4,4,4):-2}.items():
    for p in set([(i,j,k),(i,k,j),(j,i,k),(j,k,i),(k,i,j),(k,j,i)]): kap[p]=val
H=np.array([[0,0,0,-1,1],[0,0,1,1,-1],[0,1,0,0,0],[0,0,-1,0,1],[0,-1,0,1,0],[1,0,1,0,-2],
            [0,0,0,0,1],[1,0,0,0,0],[1,0,0,-2,0],[0,1,0,-1,1],[1,1,0,-2,0],[0,0,1,0,0],
            [1,0,0,-1,-1],[0,0,0,1,0]],float)
V_of=lambda t: np.einsum('ijk,i,j,k',kap,t,t,t)/6.
tau_of=lambda t: np.einsum('ijk,j,k',kap,t,t)/2.
def Gtau(t):
    V=V_of(t); M=np.einsum('ijk,k',kap,t)
    A=-np.linalg.inv(M)/V+np.outer(t,t)/(2*V**2); return 0.5*(A+A.T),V
def flats(G,fixes,tol=1e-8):
    wq,U=np.linalg.eigh(G); S=U@np.diag(np.sqrt(wq))@U.T; Si=U@np.diag(1/np.sqrt(wq))@U.T
    A=np.array([S@f for f in fixes]); Q,_=np.linalg.qr(A.T)
    uu,ss,_=np.linalg.svd(np.eye(5)-Q@Q.T)
    B=[Si@uu[:,i] for i in range(5) if ss[i]>tol]
    return [b/np.sqrt(b@G@b) for b in B]
from scipy.optimize import fsolve
def point(Vt,rs,tw,t2,t4):
    d=-np.sqrt(rs)*Vt**(1/3.)                       # BRANCHE d = t1 - t3 < 0 (essentiel)
    def eqs(x):
        t=np.array([x[0],x[1]+d,t2,x[1],t4]); return [V_of(t)/Vt-1.,tau_of(t)[3]/(tw*Vt**(2/3.))-1.]
    x=fsolve(eqs,[176.*(Vt/1e6)**(1/3),54.34*(Vt/1e6)**(1/3)],xtol=1e-14)
    t=np.array([x[0],x[1]+d,t2,x[1],t4])
    return t if np.max(np.abs(eqs(x)))<1e-9 else None
wg=1.6956; tt=np.linspace(0,2*np.pi,400000)
def projections(t):
    tau=tau_of(t); G,V=Gtau(t); Gi=np.linalg.inv(G)
    B=flats(G,[t/(2*V),np.eye(5)[1],np.eye(5)[3]])[:2]
    out={}
    for D,lab in ((0,"K3_prime1"),(2,"dP7_plate_prime6"),(4,"Wilson_prime8"),(3,"Wilson_Dterme_prime7")):
        cov=np.eye(5)[D]/tau[D]; nr=np.sqrt(cov@Gi@cov)
        ps=[wg*(cov@b)/nr for b in B]
        cb=np.abs(np.cos(tt)*ps[0]+np.sin(tt)*ps[1])
        out[lab]=dict(p1=float(ps[0]),p2=float(ps[1]),norme=float(np.hypot(*ps)),
                      fraction_fenetre=float(np.mean((cb>=0.30)&(cb<=0.61))))
    return out,tau,G,V,len(B)
t=point(1e6,8.6411e-12,1.51,40.,62.)
pr,tau,G,V,nB=projections(t)
inv_V=[]
for Vt in (1e6,1e9,1e12,1e15,1e18):
    sc=(Vt/1e6)**(1/3.); tv=point(Vt,8.6411e-12,1.51,40.*sc,62.*sc)
    inv_V.append(dict(volume=Vt,norme_p_K3=projections(tv)[0]["K3_prime1"]["norme"]))
dep=[]
for t2,t4 in ((28.,62.),(40.,62.),(54.,62.),(40.,78.)):
    tv=point(1e6,8.6411e-12,1.51,t2,t4)
    if tv is None: continue
    p2_,_,_,_,_=projections(tv)
    dep.append(dict(t2=t2,t4=t4,dans_cone=bool(np.all(H@tv>=-1e-8)),
                    norme_p_K3=p2_["K3_prime1"]["norme"],fraction=p2_["K3_prime1"]["fraction_fenetre"]))
save("06_candidat_4686",dict(
  geometrie=dict(h11=5,h21=81,chi=-152,structure="dP7 dP7 K3 W W W",base=[1,5,6,7,8]),
  kappa_non_nuls={"(0,3,4)":2,"(1,1,1)":2,"(1,1,3)":-2,"(1,3,3)":2,"(2,2,2)":2,"(2,2,4)":-2,
                  "(2,4,4)":2,"(3,3,3)":-2,"(3,3,4)":-2,"(4,4,4)":-2},
  hyperplans_cone=H.tolist(), branche="d = t1 - t3 < 0",
  point=dict(t=t.tolist(),tau=tau.tolist(),volume=float(V),
             eig_G=np.linalg.eigvalsh(G).tolist(),definie_positive=bool(np.linalg.eigvalsh(G).min()>0),
             dans_cone=bool(np.all(H@t>=-1e-8)),rang_complement=nB),
  projections=pr, invariance_volume=inv_V, dependance_t2_t4=dep,
  tadpoles_D3=dict(sigma_7=dict(N_O3=6,Q_D3=2),sigma_8=dict(N_O3=14,Q_D3=4))))

# ============================================================ 7. POTENTIEL DE BOUCLES
print("[7] potentiel de boucles : test de generictie")
t4min,t4max=54.337,121.774; grid=[]
for t4 in np.linspace(t4min*1.01,t4max*0.99,40):
    lo=max(t4-54.337,0.4)
    for t2 in np.linspace(lo*1.01+0.2,t4*0.99,40):
        tv=point(1e6,8.6411e-12,1.51,t2,t4)
        if tv is not None and np.all(H@tv>=-1e-8): grid.append((t2,t4,tau_of(tv),tv))
TA=np.array([g[2] for g in grid]); T2=np.array([g[0] for g in grid]); T4=np.array([g[1] for g in grid])
rng=np.random.default_rng(5); nint=0; ps=[]
for _ in range(40):
    C=rng.uniform(0.1,2.,5); i=int(np.argmin((TA*(C**2)).sum(axis=1)))
    lo=max(T4[i]-54.337,0.4)
    if T4[i]>t4min*1.06 and T4[i]<t4max*0.94 and T2[i]>lo*1.10+0.5 and T2[i]<T4[i]*0.96:
        nint+=1; ps.append(projections(grid[i][3])[0]["K3_prime1"]["norme"])
ps=np.array(ps) if ps else np.array([np.nan])
save("07_boucles",dict(points_grille=len(grid),minima_interieurs=f"{nint}/40",
  norme_p_K3_aux_minima=dict(min=float(np.nanmin(ps)),max=float(np.nanmax(ps)),mediane=float(np.nanmedian(ps))),
  fraction_dans_fenetre=float(np.mean((ps>=0.30)&(ps<=0.61))),
  formes="KK: +(g_s C_i)^2 tau_i/V^2 ; winding: -C_a/(V t_a)",
  limite="les coefficients C^KK et C^W ne sont pas calculables pour une CY generale"))
print("\nTOUS LES JSON ECRITS DANS 06_preuves/")
