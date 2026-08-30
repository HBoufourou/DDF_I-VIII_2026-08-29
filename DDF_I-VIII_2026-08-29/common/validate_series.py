"""TEST DE VALIDATION v2 — série finale CSG après les 15 corrections (04/08/2026).
Trois niveaux : (A) les scripts de dossier reproduisent les nombres ;
(B) les 15 corrections sont présentes dans les PDF ; (C) cohérence inter-articles."""
import subprocess,re,sys
def pdf(p):
    t=subprocess.run(['pdftotext','-layout',p,'-'],capture_output=True,text=True).stdout
    t=re.sub(r'\s+',' ',t)
    return t.replace('\u2019',"'").replace('\u2018',"'")
P={ 'I':pdf('article_I_stage/05_pdf_final/Boufourou_2026_DDF_I_the_stage.pdf'),
    'II':pdf('article_II_well_tower/05_pdf_final/Boufourou_2026_DDF_II_well_tower.pdf'),
    'III':pdf('article_III_medium/05_pdf_final/Boufourou_2026_DDF_III_the_medium.pdf'),
    'IV':pdf('article_IV_shadow/05_pdf_final/Boufourou_2026_DDF_IV_one_scale_shadow.pdf')}
tests=[]
def T(name,cond): tests.append((name,bool(cond)))
# ---------- A. scripts -> nombres ----------
def run(cmd): return subprocess.run(cmd,shell=True,capture_output=True,text=True).stdout
a1=run("python3 article_I_stage/02_scripts/core_scale.py")
T("A1 core_scale : M5=3.6e8 (piR), especes 1.5e9, piR=25.8",("3.6e+08" in a1) and ("1.5e+09" in a1) and ("25.8" in a1))
a2=run("python3 article_II_well_tower/02_scripts/u_derive.py")
T("A2 u_derive : u=4.434, T=9.8 TeV, w1=5.48, cb=0.30",all(k in a2 for k in["4.434","9.8 TeV","5.48","0.30"]))
a3=run("python3 article_III_medium/02_scripts/medium_numbers.py")
T("A3 medium : Bose 7.7e12, Soleil 18 kau, 1.9 TeV",all(k in a3 for k in["7.7e+12","18 kau","1.9 TeV"]))
a5=run("python3 article_I_stage/02_scripts/radion_relic.py")
T("A5 radion_relic : T_osc=4.2 TeV, dR/R<7.4e-08, 3.6e+12",all(k in a5 for k in["4.2 TeV","7.4e-08","3.6e+12"]))
a4=run("python3 article_IV_shadow/02_scripts/shadow.py")
T("A4 shadow : mu=24.1, accord 29%, 0.093",all(k in a4 for k in["24.1","29%","0.093"]))
# ---------- B. les 15 corrections dans les PDF ----------
B=[("I P5-1 mécanisme requis déclaré",'I',"which this framework does not yet contain"),
("I P5-1b estimation-pas-prédiction",'I',"an order estimate, not a prediction"),
("I P5-1c plancher données (partie 1)",'I',"the data set the floor"),
("I P5-1c fenêtre théorie (partie 2)",'I',"theory sets the window"),
("I P5-2 paragraphe fenêtre→raie",'I',"What would convert the window into a line"),
("I P5-3 exécuteur (v)",'I',"(v) The radion window"),
("II P1-1 « most severe »",'II',"the framework's most severe naturalness tension"),
("II P1-1b « locates it »",'II',"does not resolve this tension; it locates it"),
("II P1-2 interrogatoire enregistré",'II',"The interrogation must be recorded"),
("III P4-1 structure précise (2 hypothèses)",'III',"cannot distinguish two hypotheses"),
("III P4-2 trois routes",'III',"Three routes are available"),
("III P4-3 exécuteur (iv) double sens",'III',"both outcomes are scheduled by the data"),
("IV P3-1 « not derivation of the relation »",'IV',"not derivation of the relation"),
("IV P3-2 objection anticipée",'IV',"An honest objection must be met here"),
("IV P3-3 convergence = évidence",'IV',"it is the evidence for it"),
("IV P2-1 colonne Provenance",'IV',"Provenance"),
("IV P2-2 « guest in the spectrum »",'IV',"a guest in the spectrum, not a resident"),
("IV P2-3 gravitino hors convergence",'IV',"plays no role in it"),
("IV P1-1 « launder its declared debts »",'IV',"launder its declared debts"),
("IV P1-2 exécuteur (9)",'IV',"(9) The source tadpole"),
("I v2.1 troisième verrou (§5.3)",'I',"The third lock"),
("I v2.1 borne 1 partie sur 1e7",'I',"one part in ten million"),
("I v2.1 fine-tuning déclaré",'I',"filed [Open, declared] rather than dismissed"),
("I v2.1 remarque structurelle → IV",'I',"the scale of the wall that cast it"),
("I v2.1 exécuteur (vi)",'I',"(vi) The relic")]
for n,art,k in B: T("B "+n, k in P[art])
# ---------- C. cohérence inter-articles ----------
T("C1 M5=3.6e8 (I, piR) & jamais 2.8e8 nulle part",("3.6" in P['I']) and not any("2.8×108" in P[x].replace(" ","") for x in P))
T("C1b échelle des espèces dans I (raccord littérature)",("species-scale convention" in P['I']) and ("1.5" in P['I']) and ("GeV band" in P['I']))
T("C2 9.8 TeV : II (route grav) ET IV (convergence, thirty percent)",("9.8 TeV" in P['II']) and ("9.8" in P['IV']) and ("thirty percent" in P['IV']))
T("C3 7.6 TeV : IV",("7.6" in P['IV']))
T("C4 w1 : II (5.5) cohérent avec I (formule)",("5.5" in P['II']) and ("2w1" in P['I'].replace(" ","") or "w1" in P['I']))
T("C5 4.3 meV : II (cadran) ET IV (spectre)",("4.3" in P['II']) and ("4.3" in P['IV']))
T("C6 tension g~2e3√M5 : II ET IV (renouvelée)",("2×103" in P['II'].replace(" ","")) and ("2×103" in P['IV'].replace(" ","")))
T("C7 fenêtre 5-40 meV : I ET IV",("5–40" in P['I'] or "5--40" in P['I']) and ("5–40" in P['IV'] or "5--40" in P['IV']))
T("C8 51/41 ordres : II ET IV",all(("fifty-one" in P[x]) and ("forty-one" in P[x]) for x in ('II','IV')))
T("C9 0.69 (dégénérescence) : III",("0.69" in P['III']))
T("C10 une décade falsifiable : IV",("one decade" in P['IV']))
T("C11 T_osc ~ M_S : I (4 TeV) cohérent avec la bande IV (4-10)",("4 TeV" in P['I']) and ("4–10" in P['IV'] or "4-10" in P['IV']))
# ---------- rapport ----------
npass=sum(1 for _,c in tests if c)
print(f"{'TEST':64s} RESULTAT")
print("-"*76)
for n,c in tests: print(f"{n:64s} {'PASS' if c else '### FAIL ###'}")
print("-"*76)
print(f"BILAN : {npass}/{len(tests)} PASS")
sys.exit(0 if npass==len(tests) else 1)
