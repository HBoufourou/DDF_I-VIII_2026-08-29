# INFRASTRUCTURE COMMUNE — livraison 1/3

Trois fichiers partagés par les sept articles. **Testés ensemble** : compilation sans
erreur, bibtex sans erreur, les macros, le bloc de série, l'appendice et la déclaration
d'outils rendent correctement.

## `ddf-common.tex` — le préambule

À mettre en tête de chaque article par `\input{ddf-common}`. Il porte :

**La numérotation de série**, en une macro. `\seriesdate{III}{The medium}` produit
*Article III of VII · The medium · <date>*. **Change « VII » dans ce seul fichier et les
sept articles suivent** — plus jamais de « of V » et « of VI » qui traînent.

**Les étiquettes épistémiques** : `\Derived`, `\Input`, `\Posited`, `\Candidate`,
`\Open`, `\Verified`. Uniformes dans les sept.

**`\toolsstatement`** — la déclaration d'usage d'outils, identique partout.

**`\traceability{V}{note}`** — le bloc de traçabilité, avec le lien GitHub et la
structure de dossier.

## `ddf-refs.bib` — 34 entrées, clés stables

Les sept articles sont `DDF1`..`DDF7`, avec leurs titres définitifs. Les externes
gardent les clés que tes articles utilisent déjà — `MVV`, `CPV`, `Shukla`, `CCQ`, `GM`,
`BHP`, `CBQ` — plus celles qui manquaient : `SPARC`, `McGaugh16`, `Clowe`, `vanDokkum`,
`Lee20`, `EotWash`, `CMBS4`, `Planck18`, `BK18`, `CQSV`.

**Une seule bibliographie pour tout le corpus.** Une référence corrigée une fois l'est
partout.

## `appendix-where.tex` — l'appendice de liaison

C'est lui qui répond à ta question, et c'est le morceau qui relie les quatre premiers
articles aux trois derniers. Cinq sections :

**A.1 — les trois lieux.** Un tableau : qui vit sur la brane, dans le bulk, sur la paroi
lointaine, et dans quel article c'est établi.

**A.2 — la matière noire.** Le réservoir de la tour (bulk, Articles I-IV) et la matière
projetée (brane, Article VI) sont **deux énoncés différents, et l'appendice le dit**.
Ce qu'ils partagent — froid, sans collision, portant la masse — suffit à la
phénoménologie ; ce qui diffère est le constituant. **Ils sont compatibles si le
réservoir n'est pas requis d'être toute la matière noire**, et l'argument de l'Article III
demande seulement qu'il domine le budget. Étiqueté `[Open, declared]`.

**A.3 — l'énergie sombre.** L'échec du moteur scalaire à quarante ordres est **enregistré,
pas caché**, et la relocalisation vers la courbure extrinsèque est déclarée avec son
réglage.

**A.4 — la RAR et ce qui fait tourner les galaxies.** La section que tu demandais. Elle
dit la chose en une phrase :

> **Les galaxies tournent à cause du médium, elles pèsent à cause du réservoir.**

Les baryons sont sur la brane, ils sourcent Φ qui se propage dans le bulk avec une portée
de 8,2 μm — microscopique, donc aucune force entre corps de laboratoire. Ce qui agit à
l'échelle galactique n'est pas le champ libre mais le **milieu** : la composante cohérente
répond à la distribution baryonique, et cette réponse, relue sur la brane, est
l'accélération supplémentaire. Elle **sature**, d'où une relation avec une échelle g†
plutôt qu'un rescaling de G. Elle couple à la trace, donc à la masse et non à la
composition. Elle est portée par ≤ 1,4 % du secteur, donc **elle ne concurrence jamais le
réservoir dans le budget de masse**. Et elle **meurt** en champ faible (les UDG) et dans
la bulle de décohérence solaire (Cassini satisfaite identiquement).

**Et surtout : aucune des deux relocalisations ne la touche.** La réponse galactique n'a
jamais été assignée au secteur d'énergie sombre, et ce n'est pas le réservoir non plus —
c'est le troisième rôle, celui que le cadre a été construit pour expliquer.

**A.5** — ce qu'un lecteur doit emporter d'un article à l'autre : I-IV tiennent seuls,
V-VII ajoutent une origine UV, et les deux relocalisations ne sont pas rétroactives.

## Ce qui suit

**Livraison 2** : les dossiers `article-1/` à `article-4/` — sources migrées sur le
préambule commun, scripts, JSON, CSV, figures, PDF, VERIFICATION.txt.
**Livraison 3** : `article-5/` à `article-7/`, plus le README GitHub racine et le
paquet Overleaf.
