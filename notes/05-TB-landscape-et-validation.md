# TB pour le TBG — paysage physique et chaîne de validation

Note de passation rédigée avant le brainstorm sur la forme du modèle tight-binding.
Objectif macro : un TB sur une vraie cellule de graphène bicouche twisté (TBG),
alimenté par des Wanniers monocouche compressés sur gaussiennes.

## 1. Les trois familles de modèles pour le TBG

**(1) Modèles continus (Bistritzer–MacDonald, BM).** Standard bas-énergie,
c'est ce que la thèse dérive ([CGG23], [BM11]). Modèle de Dirac k·p : chaque
couche apporte un cône, couplés par un tunneling moiré-périodique (amplitudes
w_AA, w_AB). Angle magique (~1.1°) quand α = w/(v_F k_θ) ≈ 0.586. **Modèle à
un corps**, prédit déjà les bandes plates. Ce n'est PAS un TB atomistique.

**(2) TB atomistique (Slater–Koster).** Un p_z par atome, hopping
distance/orientation-dépendant entre paires. Cellule moiré angle magique
≈ 11 000 atomes → Hamiltonien creux ~10⁴×10⁴ par point k. **C'est là que la
compression gaussienne sert** : remplacer les hoppings empiriques par des
⟨w_i|H|w_j⟩ calculés vite via gaussiennes.

**(3) TB par Wannier (notre approche).** Piège connu : on ne peut PAS
wannieriser directement les deux bandes plates du TBG (obstruction de *fragile
topology* aux points de Dirac → modèles à 5/8/10 bandes chez les autres).
**Notre approche contourne ça** : on part du Wannier p_z **monocouche** (bien
défini, sans obstruction), on l'approxime par gaussiennes, et la base bicouche
≈ union des Wanniers des deux couches. On ne wannierise jamais les bandes
plates directement.

## 2. Décision retenue : un corps d'abord

**Les bandes plates de l'angle magique sont un phénomène à un corps** (BM,
non-interagissant, les prédit). Donc :

| Niveau | Capture | Intégrales requises |
|---|---|---|
| **Un corps** `H = -½∇² + V_KS` | structure de bandes, **bandes plates**, minibands moiré | overlap ✓, kinetic (à faire, facile), V_KS (déjà Fourier) |
| **Deux corps** (+Hubbard U) | isolants corrélés, supraconductivité | ERI bi-électroniques (Boys + McMurchie-Davidson) |

**Décision (validée avec l'utilisateur)** : on vise le **un corps** pour le
résultat phare (bandes plates). Le code actuel
(`src/tight_binding/hamiltonian_scalar_product.jl`, termes
`[:kinetic, :PBE_local_potential]`) est déjà sur ce rail. Les termes deux corps
sont un **plus pour une future contribution, à garder pour bien plus tard** —
hors scope immédiat.

## 3. Conséquence sur les intégrales (rappel phase « GaIn-free »)

La chaîne TB un corps n'a presque pas besoin de réimplémenter GaIn :

| Intégrale | Statut | Backend |
|---|---|---|
| `:overlap` | ✅ natif (phase A) | `src/integrals/julia_integrals.jl` |
| `:kinetic` (= ½ `:laplacian`) | ❌ à faire | actuellement GaIn |
| `:PBE_local_potential` | ✅ déjà GaIn-free | Fourier (`potential_scalar_prod`) |
| `:coulomb`, `:atomic` (ERI) | non utilisés par la chaîne TB | (deux corps, plus tard) |

**Seule pièce native manquante : `laplacian_julia`.** Facile car le Laplacien
d'une gaussienne-polynôme reste une gaussienne-polynôme :

  ∇²(p·g) = [∇²p − 4ζ(r−R)·∇p + (4ζ²|r−R|² − 6ζ)·p] · g  =: q(r)·g

donc ⟨X₁, ∇²X₂⟩ = recouvrement entre X₁ et la nouvelle gaussienne-polynôme q·g,
que `overlap_julia` sait déjà calculer. Pas de fonction de Boys, pas de
McMurchie-Davidson. Estimation ~½ journée avec tests. Auto-validable contre la
grille (appliquer ½|G|² en Fourier sur X₂, dotter avec X₁ — machinerie existante,
cf. `_check_L²_dot_precision`).

## 4. Échelle de validation (du plus simple au plus probant)

Chaque marche valide la précédente :

1. **Monocouche** : retrouver les cônes de Dirac, la vitesse de Fermi correcte
   (sanity du Wannier p_z unique).
2. **Bicouche non-twistée** (AB/AA) : bandes bicouches simples — vérifie le
   **hopping interlayer**.
3. **TBG commensurable grand angle** (petite cellule) : minibands moiré.
4. **Près de l'angle magique** : aplatissement des bandes. *Le* résultat.

À chaque marche, la validation « chaîne complète » = comparer ⟨w_i|H|w_j⟩
calculé avec le **vrai Wannier plane-wave** vs le **Wannier gaussien-compressé**.
L'erreur sur les hoppings → erreur sur les bandes.

## 5. Drapeau pratique : la queue p_z hors-plan

Le hopping **interlayer** ⟨w_i^haut | H | w_j^bas⟩ est l'ingrédient moiré crucial
et dépend de la **queue p_z en z** (séparation ~3.35 Å). Si la compression
gaussienne reproduit mal cette queue, les bandes plates seront fausses même si
le fit H¹ global est bon.

→ **Hypothèse à challenger au brainstorm** : la métrique de validation
pertinente n'est peut-être pas l'erreur H¹ globale (phase B) mais **l'erreur
directe sur les hoppings interlayer**. Possible reformulation de l'objectif de
compression orientée « hopping » plutôt que « norme H¹ ».

## 6. Questions ouvertes pour le brainstorm

1. Forme exacte du TB un corps : quels hoppings inclure (cutoff de distance),
   on-site vs intra/inter-layer, comment échantillonner la cellule moiré ?
2. La métrique de compression doit-elle passer de « erreur H¹ » à « erreur sur
   les hoppings » (cf. §5) ? Si oui, ça boucle sur la phase B (loss à redéfinir).
3. Échelle des expériences : à quel angle commensurable commencer pour avoir une
   petite cellule tractable en local avant de viser l'angle magique ?
4. Le terme V_KS : `potential_scalar_prod` exige le `scfres` (potentiel KS DFT).
   Disponible pour la monocouche ? Comment l'obtenir pour la config bicouche ?

## 7. Tâches identifiées (à ordonnancer au brainstorm)

- **T-lap** : `laplacian_julia` (native, ~½ j, auto-validable contre grille) —
  débloque `:kinetic` analytique et l'overlap H¹ analytique `Hˢ_overlap(Ms; s=1)`.
- **T-grid** : chaîne de validation sur grille (real-Wannier vs Gaussien) sans
  aucune lib, étape 1–2 de l'échelle §4. Sert de référence pour tout analytique.
- **T-band** : reconstruction de la structure de bandes monocouche depuis le TB
  compressé, comparaison aux bandes DFT.
- (plus tard) montée en cellule TBG, puis éventuellement deux corps.

## 8. Statut du reste du projet (rappel)

- Phase A (conditionnement) : terminée, voir `02-bilan-phase-A.md`.
- Phase B (fit conjoint ML) : Tasks 1–9/12 faites, voir `RESUME-HERE.md`.
  Décision en suspens : option (C+A) run scaled-down ForwardDiff vs (B) Enzyme.
  Cette chaîne de validation TB est **un thread parallèle** à la phase B :
  on peut la mener indépendamment (elle ne dépend que du Wannier compressé,
  qu'on a déjà via le greedy de phase A).

## 9. Références

- [BM11] Bistritzer, MacDonald. *Moiré bands in twisted double-layer graphene*. PNAS 2011.
- [CGG23] Cancès, Garrigue, Gontier. *Simple derivation of moiré-scale continuous models for TBG*. PRB 2023.
- [Cao+18] Cao et al. *Unconventional superconductivity in magic-angle graphene superlattices*. Nature 2018.
- Manuscrit Vidal 2024, chap. 5 (BM, CGG, compression Wannier→GTO).
- Slater-Koster atomistique : Trambly de Laissardière et al. ; Moon, Koshino.
- Fragile topology / obstruction Wannier : Po-Zou-Vishwanath-Senthil ; Kang-Vafek.
