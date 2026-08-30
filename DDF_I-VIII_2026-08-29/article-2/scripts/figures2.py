import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal
def tour(u,mu=0.0,N=3000):
    h=u/N; x=np.linspace(h,u-h,N-1)
    V=np.sinh(mu*x)/mu if mu>0 else x
    d=1/h**2+V; e=-1/(2*h**2)*np.ones(N-2)
    Eo=eigh_tridiagonal(d,e,select='i',select_range=(0,1))[0]
    d2=d.copy(); d2[0]=1/(2*h**2)+V[0]
    Ee=eigh_tridiagonal(d2,e,select='i',select_range=(0,1))[0]
    E=np.sort(np.concatenate([Ee,Eo]))[:4]
    return E/E[0]
u0=4.434; hbar_c=197.3; ell=5.809; piR=25.76
# fig1 : profil + tour
fig,(a1,a2)=plt.subplots(1,2,figsize=(9.8,3.9))
y=np.linspace(0,u0,200)
a1.plot(y,1-y/u0,color='#1a5276',lw=2.4)
a1.annotate('$+J_0$',(0.05,1.0),fontsize=12,color='#c0392b'); a1.annotate('$-J_0$',(u0-0.35,0.03),fontsize=12,color='#c0392b')
a1.axvline(0,color='k',lw=3); a1.axvline(u0,color='grey',lw=3)
a1.set_xlabel('$y/\\ell$'); a1.set_ylabel('$\\sigma(y)$ (arb.)'); a1.set_title("(a) $\\sigma''=0$: exactly linear; Gauss forces $J_\\pi=-J_0$")
x=np.linspace(0.01,u0,300); V=x
a2.plot(x,V,'k--',lw=1)
E=tour(u0); cols=['#0e6251','#1a5276','#7d3c98','#c0392b']
a4=[2.338]  # normalisation visuelle
for i,(En,c) in enumerate(zip(E,cols)):
    Ev=En*2.338
    a2.hlines(Ev,0,min(Ev,u0),color=c,lw=2.2)
    a2.annotate(f'{En:.3f}',(min(Ev,u0)+0.06,Ev),fontsize=9,color=c,va='center')
a2.axvline(0,color='k',lw=3); a2.axvline(u0,color='grey',lw=3)
a2.set_xlim(-0.15,5.6); a2.set_xlabel('$y/\\ell$'); a2.set_ylabel('$E/\\varepsilon$'); a2.set_title(f'(b) interleaved tower, $u={u0}$')
plt.tight_layout(); plt.savefig('01_source/figures/fig1_welltower.png',dpi=170); plt.close()
# fig2 : les deux cadrans
fig,(b1,b2)=plt.subplots(1,2,figsize=(9.6,3.8))
ms=np.linspace(0.01,6,26); r0=tour(u0)
dev=[np.max(np.abs(tour(u0,m*ell/hbar_c)-r0)/r0)*100 for m in ms]
b1.plot(ms,dev,color='#1a5276',lw=2)
b1.axhline(1,color='#c0392b',ls='--',lw=1); b1.axvline(4.3,color='#c0392b',ls=':',lw=1)
b1.annotate('1% on ratios',(0.2,1.12),fontsize=9,color='#c0392b')
b1.annotate('$m_\\sigma=4.3$ meV',(4.35,0.15),fontsize=9,color='#c0392b',rotation=90)
b1.set_xlabel('$m_\\sigma$ [meV]'); b1.set_ylabel('max ratio deviation [%]'); b1.set_title('(a) dial 1: the source mass')
us=np.linspace(4.1,4.8,20); r4=[tour(x)[3] for x in us]
b2.plot(us,r4,color='#0e6251',lw=2)
b2.plot([u0],[tour(u0)[3]],'o',ms=9,color='#c0392b')
b2.annotate('derived $u=4.43$',(u0+0.02,tour(u0)[3]),fontsize=9,color='#c0392b')
b2.set_xlabel('$u=\\pi R/\\ell$'); b2.set_ylabel('fourth ratio $E_4/E_1$'); b2.set_title('(b) dial 2: the fourth ratio measures $u$')
plt.tight_layout(); plt.savefig('01_source/figures/fig2_dials.png',dpi=170)
print("figures II OK")
