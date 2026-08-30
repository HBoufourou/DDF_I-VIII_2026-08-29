# Paquet Overleaf — les sept articles dans un seul projet

Téléverse ce dossier tel quel dans un projet Overleaf vierge.

**Sept documents indépendants** : `article-1.tex` … `article-7.tex`. Chacun compile seul.
Dans Overleaf, choisis le fichier principal en haut à gauche (*Main document*).

**Trois fichiers partagés**, à ne pas dupliquer :
`ddf-common.tex` (préambule, macros, numérotation de série, déclaration d'outils,
bloc de traçabilité), `ddf-refs.bib` (une bibliographie pour les sept),
`appendix-where.tex` (l'appendice de liaison).

**Réglage Overleaf** : compilateur *pdfLaTeX*, et coche *Normal* pour que bibtex tourne.
Si la bibliographie ne s'affiche pas au premier passage, recompile — bibtex a besoin de
deux tours.

**Changer la longueur de la série** : une seule ligne dans `ddf-common.tex`, la macro
`\seriesdate`. Les sept suivent.
