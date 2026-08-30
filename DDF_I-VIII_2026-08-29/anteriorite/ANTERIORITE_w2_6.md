# Dossier d'antériorité (R11) — |w⃗|² = 6
### Article V, §3 — préalable au dépôt arXiv des Articles I–VII
**Établi le :** 23 août 2026 · **Script :** `verif_normes_w2_6.py` (toutes assertions passées)

---

## 1. La revendication telle qu'imprimée

Article V, §3 : réduction séquentielle 10→9 (intervalle) puis 9→4 (X5),
poids canoniques ×√2, vecteur tension-couplage de la DBI D8

  w⃗ = (5√2/4, √(10/7), 9/(2√14)),  |w⃗|² = 25/8 + 10/7 + 81/56 = **6**

avec |v⃗|² = 14 (terme F₀), |u⃗|² = 24 (mur), w⃗·v⃗ = 8, w⃗·u⃗ = 11, et la phrase :

> « The identity |w⃗|² = 6 holds for every Dp transverse to the first reduced
> block; an anteriority search returned nothing, and it is offered here for
> correction rather than as a priority claim. »

## 2. Journal de recherche (exécuté le 23/08/2026)

| # | Requête | Résultat utile |
|---|---|---|
| R1 | invariant Δ, couplage dilatonique des p-branes, réduction dimensionnelle | a² = Δ − 2dd̃/(D−2) ; **Δ préservé par réduction KK sans troncature** — hep-th/9605082 (« Vertical versus Diagonal Dimensional Reduction for p-branes »), hep-th/9508042 (CTP TAMU-31/95 / Imperial TP-94-95/56) ; restaté avec formules de transport dans arXiv:1603.08084, Eqs. (54)–(55). *(Auteurs exacts à confirmer au moment de la citation : lignée Lü–Pope(–Stelle).)* |
| R2 | exposants du potentiel scalaire (dilaton, volume) pour flux / D-branes / O-plans en 4D | Tableau standard à deux modules — arXiv:1003.0029 (Wrase–Zagermann) : V_H ∝ τ⁻²ρ⁻³, V_Fp ∝ τ⁻⁴ρ^(3−p), **V_Oq/Dq ∝ τ⁻³ρ^((q−6)/2)** ; lignée Hertzberg–Kachru–Taylor–Tegmark (0711.2512). |
| — | Garousi–Myers (déjà cité [8] de l'Article V) | coefficient dilatonique de disque (p−3)/2√2 : identité confirmée composante par composante. |

Aucune des sources ne contient l'énoncé « |w⃗|² = 6 » ni la base à trois
modules (dilaton, volume X5, intervalle).

## 3. Ce que le calcul établit (script, zéro paramètre)

**(a) |v⃗|² = 14 n'est pas nouveau.** Δ(F₀) = (5/2)² + 2(−1)(9)/8 = **4** en
10D ; l'invariance publiée transporte en 4D : c² = 4 + 3 = 7, soit **14** dans
la convention ×√2 du corpus. C'est l'invariant Δ = 4 de la masse de Romans,
lu en 4D. **À citer, pas à revendiquer.**

**(b) |w⃗|² = 6 échappe au théorème.** Le terme de bord (pas de facteur
e^(β₁ρ)) transporté donnerait Δ_eff = 3 − 3 = **0** — aucune valeur BPS
standard (4/N). Le « 6 » n'est donc pas dans les tables Δ, et c'est
structurel : l'invariance vaut pour les termes de volume, pas pour les
sources de bord. Cohérence interne : |v⃗ − w⃗|² = 4, le vecteur d'écart étant
exactement (facteur électrique-vs-source du dilaton, 0, √2β₁).

**(c) La décomposition qui situe la nouveauté.** Projeté sur le plan
(τ, ρ) de Wrase–Zagermann, w⃗ redonne **exactement** la norme publiée du
D8 : 9/2 + (8−6)²/6 = **31/6**. Le résidu orthogonal — la direction de
*forme* X5-contre-intervalle, absente de tout traitement à volume isotrope —
vaut **5/6**. D'où

  **|w⃗|² = 31/6 (publié) + 5/6 (l'apport de la scission de l'intervalle) = 6.**

La contribution propre de l'Article V est le second terme et l'exactitude de
la somme, pas la structure du premier.

**(d) ALERTE — la phrase « every Dp » paraît fausse.** Avec les poids de la
même convention (dilaton (p−3)/4 ; intervalle (p+1)α₁ ; géométrie
4α₂+(p−3)β₂ ; ×√2), un Dp transverse à l'intervalle enroulé sur un
(p−3)-cycle donne la forme fermée

  **|w⃗(p)|² = (p−3)²/5 − (p−3) + 6**

soit 6 ; 5,2 ; 4,8 ; 4,8 ; 5,2 ; 6 pour p = 3…8. **La valeur 6 n'est
atteinte que pour p = 8 et p = 3** — la brane qui enroule tout X5 et celle
qui n'en enroule rien. Si Hicham entendait d'autres poids pour les p
intermédiaires, ils doivent être écrits ; en l'état, la lecture naturelle de
l'article rend la phrase incorrecte, et c'est la revendication d'antériorité
elle-même qui est mal délimitée par elle.

## 4. Verdict d'antériorité

| Énoncé | Statut | Action |
|---|---|---|
| |v⃗|² = 14 | **corollaire direct du Δ = 4 publié** | citer la lignée Δ ; retirer toute nuance de nouveauté |
| contenu (dilaton, volume) de w⃗ = 31/6 | **publié** (tableau à deux modules) | citer 1003.0029 / 0711.2512 |
| scission de l'intervalle → +5/6, somme = 6 exactement | **non trouvé publié** | revendication étroite, défendable, ancêtres cités |
| « holds for every Dp » | **faux d'après le calcul (p ∈ {3, 8} seulement)** | corriger avant dépôt — validation de Hicham requise |
| coïncidence D3 ↔ D8 à 6 | non trouvée publiée | à signaler comme curiosité, sans revendication |

## 5. Texte proposé pour l'Article V, §3 (remplace la phrase citée en §1)

> Two of these numbers are inherited rather than new. |v⃗|² = 14 is the
> four-dimensional image of the Δ = 4 invariant of the Romans mass under
> untruncated Kaluza–Klein reduction [Δ-lineage refs]; and the
> (dilaton, overall-volume) content of w⃗ reproduces the standard
> two-modulus scaling V_D8 ∝ τ⁻³ρ [WZ/HKTT refs], of squared norm 31/6.
> What the present reduction adds is the third axis: separating the interval
> length from the overall volume contributes an orthogonal 5/6, and
> |w⃗|² = 31/6 + 5/6 = 6 exactly. The boundary term escapes the bulk
> Δ-invariance (its transported invariant would be Δ = 0, no BPS value),
> which is why this norm appears in no Δ table. The same counting gives
> |w⃗(p)|² = (p−3)²/5 − (p−3) + 6 for a Dp transverse to the interval
> wrapping a (p−3)-cycle: the value 6 is attained at p = 8 and p = 3 only.
> A dated search log accompanies this article; the decomposition is offered
> with its ancestors cited, for correction rather than as a priority claim.

Et §12(ii) reste valable tel quel (« if |w⃗|² = 6 is found already published,
that claim is withdrawn and the rest stands »).

## 6. Ce qui reste à faire avant que ceci débloque le dépôt

1. **Hicham valide ou réfute le comptage général en p** (§3d). S'il le
   réfute, il écrit les poids qu'il entendait, et ce dossier est recalculé.
2. Confirmer les listes d'auteurs des références Δ au moment de mettre les
   `\bibitem` (identifiants arXiv déjà consignés ci-dessus).
3. Intégrer le texte du §5 (ou sa version amendée) dans `article-5.tex`,
   joindre ce dossier et le script au dépôt Zenodo des sept articles.

*Rien dans ce dossier ne repose sur une affirmation non exécutée : chaque
nombre est une assertion du script joint.*
