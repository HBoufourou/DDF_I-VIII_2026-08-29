import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
# fig1 : la scene
fig,(a1,a2)=plt.subplots(1,2,figsize=(10,3.8))
a1.axvspan(-0.6,0,color='#E0A73C',alpha=.9); a1.axvspan(25.8,26.4,color='#7f8c8d',alpha=.55)
a1.annotate('SM brane\n(thin core ~TeV$^{-1}$)',(0.4,0.82),fontsize=9,color='#8a6d1a')
a1.annotate('far wall\n(empty)',(22.6,0.82),fontsize=9,color='#5d6d7e')
y=np.linspace(0,25.8,300)
a1.plot(y,0.55+0.25*np.cos(np.pi*y/25.8),color='#5d6d7e',lw=1.6)
a1.annotate('gravity: spreads',(9.5,0.86),fontsize=9,color='#5d6d7e')
a1.plot(y,0.30*np.exp(-y/6),color='#38CFC4',lw=2)
a1.annotate('dark scalar $\\Phi$ (Article II)',(7.5,0.24),fontsize=9,color='#1b8f87')
a1.set_xlim(-1,26.5); a1.set_ylim(0,1); a1.set_yticks([])
a1.set_xlabel('$y$ [$\\mu$m]'); a1.set_title('(a) the corridor, $\\pi R = 25.8\\ \\mu$m')
# panneau b : 13 ordres
x=np.logspace(-19.5,-4.5,400)
core=np.exp(-(x/2e-19)**2); psi=np.exp(-(x/5.8e-6)**2)
a2.semilogx(x,core,color='#E0A73C',lw=2,label='SM core $|f_0|^2$, $L_g\\sim2\\times10^{-19}$ m')
a2.semilogx(x,psi,color='#38CFC4',lw=2,label='bulk mode $|\\psi_1|^2$, $\\ell\\sim6\\ \\mu$m')
a2.legend(fontsize=8,loc='center left'); a2.set_yticks([])
a2.set_xlabel('length [m]'); a2.set_title('(b) thirteen orders: the vertex is pointlike')
plt.tight_layout(); plt.savefig('01_source/figures/fig1_stage.png',dpi=170); plt.close()
# fig2 : plan (alpha, lambda) schematique
fig,ax=plt.subplots(figsize=(7.4,4.4))
lam=np.logspace(0,2.2,200)  # um
bound=2.5e5*lam**-3.2       # frontiere schematique calibree: alpha=1 a 38.6 um
bound=bound/ (2.5e5*38.6**-3.2)  # normalise a 1 en 38.6
ax.loglog(lam,bound,color='#C4636A',lw=2,label='torsion-balance frontier (schematic, $\\alpha$=1 @ 38.6 $\\mu$m)')
ax.fill_between(lam,bound,1e4,color='#C4636A',alpha=.12)
ax.plot([8.2],[8/3],'o',ms=10,color='#1a5276',label='graviton tower: $\\alpha=8/3$ @ $8.2\\ \\mu$m  [fixed point]')
ax.axvspan(5,40,color='#8FAFD4',alpha=.3,label='radion window 5–40 $\\mu$m (floored @ 5.11 meV)')
ax.set_xlabel('range $\\lambda$ [$\\mu$m]'); ax.set_ylabel('strength $\\alpha$')
ax.set_xlim(1,160); ax.set_ylim(1e-2,1e4); ax.legend(fontsize=8,loc='upper right')
ax.set_title('the short-range landscape (schematic)')
plt.tight_layout(); plt.savefig('01_source/figures/fig2_ranges.png',dpi=170)
print("figures OK")
