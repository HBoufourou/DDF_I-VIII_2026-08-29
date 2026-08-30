# CHANGELOG — Corpus DDF, dépôt Zenodo v1.0
### 29 août 2026 · corrections appliquées aux sources avant premier dépôt

## Corrections du registre APPLIQUÉES dans cette version (sources .tex)

| Réf. | Article | Modification |
|---|---|---|
| (a) | III, §8(ii) | falsificateur corrigé : « beyond the decoherence edge » → **« within »** — le texte antérieur inversait le test Gaia DR4 |
| (b) | III (×3) + IV (écho) | mâchoire supérieure resserrée : f ≤ **8,42×10⁻³** (nul à 30 kau, arXiv:2608.24556) ; bulle solaire minimale 18 → **30 kau** ; provenance « tightened twice » (AGC 114905 puis compagnon) |
| (c) | V, §3 | « holds for every Dp » (**faux**) remplacé par le texte d'antériorité : \|v⃗\|²=14 = corollaire du Δ=4 publié ; \|w⃗\|² = 31/6 (publié) + 5/6 (scission d'intervalle) ; forme générale \|w⃗(p)\|² = (p−3)²/5 − (p−3) + 6, valeur 6 en p ∈ {3, 8} seulement |
| (d) | V | ascendances citées en ligne : hep-th/9605082, hep-th/9508042 (lignée Lü–Pope–Stelle) ; arXiv:1003.0029 (Wrase–Zagermann) — clés \cite propres à la réécriture v2 |
| (e) | VII | lacune déclarée : les sélections atteignent X produit en désintégrations de hadrons beaux ; la production **prompte** n'est pas couverte |
| (g) | IV | note de bas de page définissant le seuil mobile s_c ∝ M^{1/3} et renvoyant au protocole pré-enregistré horodaté séparément |

Le dossier d'antériorité complet (journal daté + script) accompagne le dépôt
dans `anteriorite/`, conformément au texte de l'Article V.

## NON appliqué dans v1 — assumé et daté

- **(f) Affiliation** : les articles disent « Brussels » ; à harmoniser sur
  décision de l'auteur (v1 déposée telle quelle si non tranché ce soir).
- **(j) Correction Voros (Article IX)** : bloc LaTeX prêt ; s'applique quand
  IX rejoint le dépôt (v1.1). IX ne dépend pas de LHCb.
- **Remarques T1** (jonction D⊕N comme condition de cohérence, Articles
  II/VI ; enrichissements IX) : étiquetées [Derived] mais **R11 ciblé non
  effectué** → réservées à la v2 (réécriture complète).
- **Article X** : interne, hors dépôt, par décision de programme.

## Périmètre v1

Articles I–VIII (sources, PDF, figures, données, scripts, infrastructure
commune, journal de recherche, valeurs certifiées). La v2 (réécriture
unique, registre complet) suivra les verdicts LHCb + Gaia DR4 + relectures.

## ⚠️ AVANT TÉLÉVERSEMENT — quatre PDF à régénérer

Les sources de **III, IV, V, VII** ont changé ; leurs PDF joints datent
d'avant correction. Dans Overleaf (projet existant, `ddf-refs.bib` en
place) : remplacer ces quatre .tex par ceux du paquet, recompiler,
exporter, substituer les PDF dans le dossier, puis téléverser. Quinze
minutes ; la cohérence source↔PDF d'un dépôt de priorité n'est pas
négociable.
