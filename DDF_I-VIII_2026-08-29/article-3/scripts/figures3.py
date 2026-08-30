import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,(a1,a2)=plt.subplots(1,2,figsize=(9.8,3.8))
f=[0.27,0.39,0.22,0.12]
a1.bar(range(1,5),f,color=['#38CFC4','#1a5276','#1a5276','#1a5276'])
a1.annotate('condensate\n(mediator)',(0.72,0.05),fontsize=8,color='#0e6251')
a1.annotate('reservoir: 73% at birth',(1.8,0.36),fontsize=9,color='#1a5276')
a1.set_xticks(range(1,5)); a1.set_xlabel('tower level'); a1.set_ylabel('initial fraction  [Posited]')
a1.set_title('(a) occupation at birth (flat displacement)')
E=[1,2.294,3.189,4.029]
for i,e in enumerate(E): a2.hlines(e,0.2,0.8,color='#1a5276',lw=2.5)
for i in range(3,0,-1):
    a2.annotate('',xy=(0.5,E[i-1]+0.06),xytext=(0.5,E[i]-0.06),arrowprops=dict(arrowstyle='->',color='#c0392b',lw=1.6))
a2.annotate('scalar cascade\n(Bose-stimulated $\\times 10^{10-12}$)',(0.56,2.4),fontsize=8,color='#c0392b')
a2.annotate('graviton decay: $e^{-10^{24}}$ (dead)',(0.1,4.35),fontsize=8,color='#7f8c8d')
a2.hlines(E[0],0.2,0.8,color='#38CFC4',lw=4)
a2.annotate('condensate $\\to$ vertex $\\to$ baryons',(0.12,0.6),fontsize=8,color='#0e6251')
a2.set_xlim(0,1.05); a2.set_ylim(0.3,4.7); a2.set_xticks([]); a2.set_ylabel('$E/\\varepsilon$')
a2.set_title('(b) the drip-feed')
plt.tight_layout(); plt.savefig('01_source/figures/fig1_medium.png',dpi=170); plt.close()
# fig2
fig,(b1,b2)=plt.subplots(1,2,figsize=(9.8,3.8))
rho=np.logspace(-2,1.2,100); K=4.9e-4
for ff,c,lab in ((1.2e-3,'#1a5276','$f=1.2\\times10^{-3}$'),(1.4e-2,'#0e6251','$f=1.4\\times10^{-2}$')):
    b1.loglog(rho,K/(ff*rho),color=c,lw=2,label=lab)
b1.fill_between(rho,K/(1.4e-2*rho),K/(1.2e-3*rho),color='#8FAFD4',alpha=.25)
b1.axhline(10/206.265e0/1e3*206265/206265,color='k') # dummy removed below
b1.clear()
for ff,c,lab in ((1.2e-3,'#1a5276','$f=1.2\\times10^{-3}$'),(1.4e-2,'#0e6251','$f=1.4\\times10^{-2}$')):
    b1.loglog(rho,K/(ff*rho),color=c,lw=2,label=lab)
b1.fill_between(rho,K/(1.4e-2*rho),K/(1.2e-3*rho),color='#8FAFD4',alpha=.25)
b1.axhline(10*1e3/206265,color='#c0392b',ls='--',lw=1)
b1.annotate('wide-binary switch: 10 kau',(0.013,10*1e3/206265*1.15),fontsize=8,color='#c0392b')
b1.plot([0.4],[K/(1.4e-2*0.4)],'*',ms=13,color='#E0A73C')
b1.annotate('Sun, local $\\rho$',(0.45,K/(1.4e-2*0.4)*0.75),fontsize=8,color='#8a6d1a')
b1.set_xlabel('environment density $\\rho$ [GeV cm$^{-3}$]'); b1.set_ylabel('$R_{\\rm dec}(1\\,M_\\odot)$ [pc]')
b1.legend(fontsize=8); b1.set_title('(a) the off-switch, $R_{\\rm dec}=K M^{1/3}/(f\\rho)$')
ups=[0.3,0.5,0.69,0.7,0.9]; daa=[-18.5,-9.3,0.0,0.8,11.7]
b2.plot(ups,daa,'o-',color='#1a5276',lw=1.8,ms=6)
b2.axhline(0,color='k',lw=0.8); b2.axvspan(0.59,0.78,color='#8FAFD4',alpha=.3)
b2.axvline(0.5,color='#c0392b',ls=':',lw=1.2)
b2.annotate('standard $\\Upsilon_\\star=0.5$',(0.505,-16),fontsize=8,color='#c0392b',rotation=90)
b2.annotate('$|\\Delta\\alpha_b/\\alpha_b|<5\\%$',(0.60,8.5),fontsize=8,color='#34495e')
b2.set_xlabel('imposed $\\Upsilon_\\star$'); b2.set_ylabel('apparent $\\Delta\\alpha_b/\\alpha_b$ [%]')
b2.set_title('(b) the declared degeneracy (zero at $\\Upsilon_\\star=0.69$)')
plt.tight_layout(); plt.savefig('01_source/figures/fig2_switch.png',dpi=170)
print("figures III OK")
