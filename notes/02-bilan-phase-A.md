# Phase A — bilan

Conclusions de la phase A (diagnostic du conditionnement de la matrice de
recouvrement S au cours du greedy de compression). Voir aussi
`01-plan-A-conditionnement.md` pour les hypothèses initiales.

## TL;DR

- Tikhonov sur le solve intérieur **suffit à débloquer le greedy** (on dépasse
  la borne de la thèse : 6.8% H¹ à 15 fonctions vs 12% à 10 antérieurement).
- Mais le diagnostic révèle **deux pathologies distinctes**, dont une que je
  n'avais pas anticipée :
  1. **Escape numérique en ζ → ∞** (artefact discrétisation, pas un vrai
     problème de redondance)
  2. **Redondance diffuse** (cosine max ~20%) entre fonctions à ζ ~ 1 sur les
     axes π-bond — la vraie limite structurelle de la méthode greedy

## Setup

- Ecut = 15 Ha, kgrid = [5,5,1], d = 10 Å
- Wannier précomputé en JSON, basis_supercell reconstruit
- max_iter = 15, tol = 1e-3 (volontairement strict pour pousser au-delà)
- max_xy_order = 3, max_z_order = 3 (n_SAGTOs = 6 par itération)
- regularization = `:tikhonov`, ε_rel = 1e-8

## Évolution observée

| iter | n_funcs | err H¹ | cond(S) | σ_min  | σ_max |
|------|---------|--------|---------|--------|-------|
| 1    | 1       | 34.0%  | 1.0     | 4.30   | 4.30  |
| 4    | 4       | 23.8%  | 24.7    | 2.21   | 54.6  |
| 8    | 8       | 10.1%  | 25.1    | 2.19   | 55.0  |
| 11   | 11      | 8.6%   | 86.7    | 0.63   | 55.0  |
| **12** | **12** | **7.8%** | **2.0e7** | **2.7e-6** | **55.0** |
| 15   | 15      | 7.4%   | 4.8e7   | 1.2e-6 | 57.4  |

L'effondrement à iter 12 est brutal (5 ordres de grandeur sur σ_min) et
discontinu — pas une dégradation graduelle.

## Pathologie 1 — escape numérique en ζ

Spreads des fonctions ajoutées (pour les itérations sensibles) :

| iter | ζ        | α (centre primaire)        |
|------|----------|----------------------------|
| 11   | 1.71e-01 | origine                    |
| **12** | **5.28e+03** | π-bond proche (α ≈ 0.5) |
| 13   | 4.17e-01 | origine                    |
| 14   | 1.15e+00 | π-bond                     |
| **15** | **6.30e+02** | origine                |

Aux itérations 12 et 15, l'optimiseur extérieur (CG sur log ζ) envoie ζ à
des valeurs 100×–1000× plus grandes que celles des itérations saines
(typiquement 0.05–2.5).

**Mécanisme** :
- Une SAGTO de spread ζ = 5283 a une largeur effective $1/\sqrt{\zeta} \approx
  0.014$ Å — 100× plus petite que la liaison C-C
- Son contenu Fourier est concentré autour de $|q|^2 \sim 4\zeta \approx 20000$,
  bien au-delà de la coupure Ecut = 15 ($|q|^2 \le 30$ environ)
- Sa **représentation plane-wave est quasi-nulle**
- Mais la normalisation L² (calculée analytiquement en Julia pur) continue à
  voir la fonction comme normée à 1
- D'où la matrice S avec un élément diagonal ≈ 4e-6 (`⟨Φ_12, Φ_12⟩_{H^1}`),
  qui devient instantanément σ_min

Bref, c'est un artefact d'incohérence entre **normalisation analytique
continue** et **produit scalaire discret** ; le greedy l'exploite pour
"trouver" des optima qui ne sont pas physiquement représentables. Ce n'est
pas une vraie redondance.

**Implications** :
- Borner ζ par le haut (ζ_max ~ Ecut/4) est non négociable pour la phase B
- Augmenter Ecut décale la borne mais ne supprime pas le mécanisme
- Le seul vrai correctif structurel est d'utiliser le produit scalaire
  *discret* aussi pour la normalisation (cohérence interne)

## Pathologie 2 — redondance diffuse

Aux itérations 1–11 (régime "sain"), cond(S) croît de 1.0 à 86.7. Pas
catastrophique, mais notable : c'est la signature de la *vraie* limite
structurelle du greedy.

Cosines normalisés des paires les plus recouvrantes (sur les 15 fonctions
finales) :

| paire     | cosine | ζ_i  | ζ_j  | direction    |
|-----------|--------|------|------|--------------|
| Φ_4 ↔ Φ_14 | 0.20   | 0.84 | 1.15 | π-bond       |
| Φ_8 ↔ Φ_14 | 0.21   | 2.32 | 1.15 | π-bond       |
| Φ_2 ↔ Φ_4  | 0.22   | 0.36 | 0.84 | π-bond       |
| Φ_4 ↔ Φ_7  | 0.17   | 0.84 | 0.60 | π-bond/centre |

**Aucune paire au-delà de ~22% de cosine.** Donc pas de quasi-colinéarité
catastrophique — H1 (mon hypothèse initiale) est partiellement validée mais
surévaluée. La redondance est **collective et modérée** : plusieurs SAGTOs à
ζ ~ 0.5–2.5 le long des axes π-bond se recouvrent un peu, sans qu'aucune
paire ne sature.

C'est exactement le genre de pathologie qu'un greedy ne peut pas résoudre :
chaque ajout est localement optimal mais l'ensemble accumule de la
redondance.

## Implications pour la phase B

Trois contraintes de design directement issues du diagnostic :

1. **Borner ζ par le haut** dans la paramétrisation (ζ_max calé sur Ecut),
   pour interdire l'escape numérique.
2. **Cohérence des produits scalaires** : ou bien tout discret, ou bien tout
   continu. Le mélange actuel est responsable de la pathologie 1.
3. **Régularisation des coefficients** dans la loss : pénalité L2 `μ ‖c‖²`
   ou contrainte explicite, pour absorber la redondance diffuse sans qu'elle
   ne soit compensée par des coefficients géants.

Indirectement, le diagnostic motive aussi :

4. Considérer **Ecut = 30 ou 50** pour avoir des produits scalaires plus
   précis (les JSONs existent déjà).
5. Le greedy a une borne intrinsèque dans ce setup : **~11 fonctions /
   8.6% H¹**. Au-delà, on entre soit dans la pathologie 1 (escape), soit
   dans la pathologie 2 (redondance accumulée). Le fit conjoint de phase B
   doit pouvoir faire mieux que ça pour justifier le travail.

## Livrables phase A

- `src/integrals/julia_integrals.jl` — overlap pur Julia (commit précédent)
- `src/basis_functions/BasisFunctions.jl` — `_make_inner_solver` avec
  4 stratégies (`:none`, `:tikhonov`, `:svd_truncation`, `:pivoted_cholesky`)
- `src/compression/compression.jl` — kwargs `regularization` et
  `regularization_param` propagés
- `workflow/diagnose_conditioning.jl` — script instrumenté avec callback,
  inspection post-run des paires les plus recouvrantes
- `workflow/scan_basis_params.jl` — utilitaire pour retrouver les params
  basis ayant produit un JSON donné
- `workflow/diag_outputs/` — artefacts (CSV log + matrice de recouvrement)
