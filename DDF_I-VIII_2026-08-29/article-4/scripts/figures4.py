import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
MPl=2.4e18
fig,(a1,a2)=plt.subplots(1,2,figsize=(10,3.9),gridspec_kw={'width_ratios':[1.35,1]})
rows=[("wall tension $T^{1/4}$",2,13,8.2),("localisation $L_g^{-1}$",1,10,None),
      ("radion $V_{\\rm stab}^{1/4}$",3.5,10,None),("nucleation",5,7,6),
      ("SUSY band",4,10,None),("meV knot lifted",7.2,8.0,7.6)]
for i,(lab,lo,hi,c) in enumerate(rows):
    a1.barh(i,hi-lo,left=lo,color='#8FAFD4',alpha=.55,height=.55)
    if c: a1.plot([c],[i],'o',color='#1a5276',ms=7)
a1.axvspan(3.5,9.8,color='#E0A73C',alpha=.18)
a1.annotate('$M_S$ from the meV sector\n(3.5–9.8 TeV)',(3.7,5.55),fontsize=8,color='#8a6d1a')
a1.set_yticks(range(len(rows))); a1.set_yticklabels([r[0] for r in rows],fontsize=8)
a1.set_xscale('log'); a1.set_xlim(0.8,20); a1.set_xlabel('TeV'); a1.set_title('(a) six quantities, one band')
mu=24.1
sp=[("$m_\\Phi$ (anchor)",24,'#1a5276'),("$m_{3/2}$ (imported)",13.9,'#7d3c98'),
    ("$m_r$ window",None,'#0e6251'),("$m_\\sigma$ bound",4.3,'#c0392b'),("$\\Lambda^{1/4}$",2.25,'#7f8c8d')]
a2.axhspan(mu/np.sqrt(10),mu*np.sqrt(10),color='#8FAFD4',alpha=.25)
a2.annotate('one decade around $\\mu$',(0.04,mu*2.4),fontsize=8,color='#34495e')
a2.fill_between([0.25,0.75],5,40,color='#0e6251',alpha=.3)
for i,(lab,m,c) in enumerate(sp):
    if m: a2.hlines(m,0.15,0.85,color=c,lw=2.4)
    a2.annotate(lab,(0.88,m if m else 14),fontsize=8,color=c,va='center')
a2.hlines(mu,0.1,0.9,color='k',lw=1,ls='--'); a2.annotate('$\\mu=24.1$',(0.02,mu*1.05),fontsize=8)
a2.set_yscale('log'); a2.set_ylim(1,120); a2.set_xticks([]); a2.set_ylabel('meV')
a2.set_title('(b) the spectrum of the shadow')
plt.tight_layout(); plt.savefig('01_source/figures/fig1_convergence.png',dpi=170); plt.close()
# fig2 : spectre (repris compact) + triangle
fig,(b1,b2)=plt.subplots(1,2,figsize=(9.6,3.8))
cs=[("$m_\\Phi$",1.0),("$m_{3/2}$",0.58),("$m_r$",None),("$m_\\sigma$",0.18)]
b1.axhspan(10**-0.5,10**0.5,color='#8FAFD4',alpha=.3)
b1.fill_between([1.8,2.2],0.21,1.7,color='#0e6251',alpha=.4)
xs=[0,1,2,3]
for x,(lab,c) in zip(xs,cs):
    if c: b1.plot([x],[c],'o',ms=9,color='#1a5276')
    b1.annotate(lab,(x,2.6),fontsize=10,ha='center')
b1.annotate('$<0.18$',(3,0.13),fontsize=8,ha='center',color='#c0392b')
b1.set_yscale('log'); b1.set_ylim(0.05,4); b1.set_xticks([]); b1.set_ylabel('$c_i=m_i/\\mu$')
b1.set_title('(a) order-one coefficients, one decade')
tri_x=[0,1,0.5,0]; tri_y=[0,0,0.85,0]
b2.plot(tri_x,tri_y,color='#34495e',lw=1.6)
b2.annotate('torsion balance\n$m_r \\Rightarrow M_S=\\sqrt{m_r M_{\\rm Pl}}$',(-0.06,-0.16),fontsize=8,ha='left')
b2.annotate('collider\n$M_S$ direct (4–10 TeV)',(1.02,-0.16),fontsize=8,ha='right')
b2.annotate('LISA\n$\\beta/H\\in[10,100]$ (cond.)',(0.5,0.9),fontsize=8,ha='center')
b2.annotate('$M_S$',(0.5,0.3),fontsize=13,ha='center',color='#c0392b')
b2.set_xlim(-0.15,1.15); b2.set_ylim(-0.35,1.1); b2.axis('off'); b2.set_title('(b) the experimental triangle')
plt.tight_layout(); plt.savefig('01_source/figures/fig2_spectrum_triangle.png',dpi=170)
print("figures IV OK")
