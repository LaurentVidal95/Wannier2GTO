# Validation TB on-site — bilan

Bilan des sessions consacrées au premier maillon de la chaîne de validation
tight-binding : comparer les intégrales TB calculées avec le vrai Wannier
plane-wave et avec le Wannier gaussien-compressé.
Spec : `06-design-validation-TB-onsite.md`. Plan : `07-plan-validation-TB-onsite.md`.

## TL;DR

- `laplacian_julia` / `kinetic_julia` sont **natifs** : la chaîne TB un corps est
  entièrement GaIn-free.
- **Un bug de convention a été découvert et corrigé** (§3) : les monômes des
  gaussiennes étaient lus en coordonnées absolues par le chemin analytique et
  relatives par tout le reste. Toutes les SAGTOs π-bond étaient intégrées faux.
  Invisible aux tests existants, qui ne testaient que des fonctions centrées.
- Après correction, chemins analytique et grille **coïncident** (1.4766869 vs
  1.4767 à Ecut=50, contre 1.4308 vs 1.4767 avant).
- **L'erreur sur le cinétique on-site est le carré de l'erreur H¹** —
  superconvergence de Galerkin, vérifiée sur deux points (§4). Avec une
  prédiction testable : ça ne doit **pas** valoir pour les hoppings.

## 1. Ce qui a été livré

| Fichier | Rôle |
|---|---|
| `src/integrals/julia_integrals.jl` | `two_center_moment_1d`, `overlap_julia`, `_axis_laplacian_terms`, `laplacian_julia`, `kinetic_julia` |
| `src/integrals/integrals.jl` | `symb_to_integral` : `:laplacian`/`:kinetic` → natif |
| `src/basis_functions/GaussianPolynomials.jl` | évaluateur réel + docstring mis en convention relative |
| `src/tight_binding/hamiltonian_scalar_product.jl` | `compare_onsite_kinetic` |
| `test/test_integrals_laplacian.jl` | forme fermée, grille, **décentré degré ≥1**, symétrie, H¹, intégration |
| `workflow/validate_tb_onsite.jl` | script de comparaison, `W2G_ECUT` paramétrable |

Commits : `69104ca` → `8988502`.

### Le principe

∇²X₂ reste une gaussienne-polynôme de **mêmes** ζ et centre, donc tout se ramène
à des appels `overlap_julia` : pas de fonction de Boys, pas de
McMurchie-Davidson, et c'est **ForwardDiff-compatible**.

Convention de signe (imposée par `Hˢ_overlap(Ms; s=1)`, `Hs_scalar_prods.jl:32`) :

    laplacian_julia(X₁,X₂) = ⟨∇X₁,∇X₂⟩ = −⟨X₁,∇²X₂⟩     (forme positive)
    kinetic_julia          = ½ laplacian_julia

**Effet de bord** : `Hˢ_overlap(Ms; s=1)` (overlap H¹ analytique) tourne
désormais sans GaIn — piste d'accélération pour la phase B.

## 2. Corrections apportées au plan initial

1. **La forme fermée (2L+3)ζ du plan était fausse.** Redérivation par les
   moments 1D :

       ⟨∇g,∇g⟩ / ‖g‖² = ζ · Σᵢ (4nᵢ−1)/(2nᵢ−1)

   (2L+3)ζ n'est valable que si **tous** les exposants valent 0 ou 1 ; z²
   mélange L=2 et L=0, d'où 13ζ/3. Bonus : cette forme exerce la branche `a≥2`
   qu'un test s-only ne touche jamais.

2. **Un couple de test du plan était dégénéré** : `xz`×`yz` s'annule par parité
   → l'assertion portait sur 0 ≈ 0. Remplacé.

## 3. Le bug de convention (découverte principale)

### Symptôme

Matrices de Gram analytique vs grille en désaccord jusqu'à **91 %**.
Cas minimal :

    s(ζ=1)@0 · px(ζ=0.5)@[1,0,0]     analytique = +0.178527     grille = −0.618436

### Diagnostic

Deux conventions coexistaient pour les monômes d'une SAGTO :

| Convention | Où |
|---|---|
| **Relative**, (r−α)ⁿ | `SAGTO_fourier_transform`, `translate`, `rotate`, `enforce_D3_symmetry`, `SAGTO_basis` |
| **Absolue**, rⁿ | `overlap_julia` (phase A), `analytic_norm`, évaluateur réel, docstring |

Elles coïncident **si et seulement si** α = 0. Or `_check_L²_dot_precision` ne
teste que des fonctions centrées en 0 : le bug était structurellement invisible.

### Pourquoi la relative est la bonne

- `SAGTO_basis` construit des polynômes *symmetry-adapted* : ils ne sont
  symétriques autour de α que s'ils sont exprimés en (r−α).
- `translate` déplace le centre en gardant le polynôme — ce n'est une
  translation **que** sous la convention relative. Idem `rotate` et
  `enforce_D3_symmetry`.
- C'est la convention standard en chimie quantique (donc de GaIn, que
  `overlap_julia` remplaçait).

### Portée

- **Fausses** : toutes les intégrales analytiques impliquant une SAGTO
  décentrée de degré ≥ 1 — c'est-à-dire les fonctions π-bond, le gros de la base.
- **Saines** : la compression elle-même, qui travaille dans la représentation
  grille (auto-cohérente en relative). **Les compressés de `data/` restent
  valides** ; c'était leur *lecture analytique* qui ne l'était pas.

### Correction

`overlap_julia` passe à l'expansion binomiale à deux centres autour du
centre-produit (`two_center_moment_1d`). En coordonnées relatives, le Laplacien
**se simplifie** — le centre disparaît, trois termes seulement :

    ∂²ᵥ[vᵃ e^{−ζv²}] = [a(a−1)v^{a−2} − 2ζ(2a+1)vᵃ + 4ζ²v^{a+2}] e^{−ζv²}

Tests de non-régression ajoutés sur des paires **décentrées de degré ≥ 1** —
le trou de couverture exact.

**Validation croisée** : analytique et grille, jusqu'ici en désaccord de 3 %,
s'accordent maintenant à 4e-5 (T = 1.4766869 vs 1.4767 à Ecut=50).

## 4. Résultats on-site (après correction)

Setup : d = 10 Å, kgrid = [5,5,1], 10 fonctions de base.

| | Ecut=15 | Ecut=50 |
|---|---|---|
| erreur H¹ (compression) | 5.59 % | 12.38 % |
| ‖w_gto‖ | 0.9989 | 0.9978 |
| T_ref (grille) | 1.4721265 | 1.5051430 |
| T_gto (analytique) | 1.4781704 | 1.4766869 |
| **erreur brute** | **0.41 %** | **1.89 %** |
| erreur normalisée | 0.63 % | 1.45 % |

Le « déficit de norme de 4 % » observé avant correction n'existait pas : c'était
le bug. La norme est bonne à 0.1 %, donc brut ≈ normalisé, et **rescaler la
compressée n'apporte rien** — c'est structurel : `project_wannier_on_basis` est
une projection H¹-**orthogonale** (⟨w−g, g⟩_H¹ = −9e-16, zéro machine), donc
déjà optimale dans son span, facteur optimal = 1.000000.

### La loi : superconvergence de Galerkin

| | (erreur H¹)² | erreur cinétique observée | rapport |
|---|---|---|---|
| Ecut=15 | 0.313 % | 0.411 % | 1.31 |
| Ecut=50 | 1.533 % | 1.891 % | 1.23 |

**L'erreur sur le cinétique est le carré de l'erreur H¹**, à 25 % près, sur deux
points très différents. Ce n'est donc pas que « H¹ est conservatif » : c'est que
pour une fonctionnelle quadratique associée au produit scalaire de la
projection, l'erreur est en ‖δ‖² et non ‖δ‖.

### Prédiction testable (importante)

Cette superconvergence tient parce que le cinétique **est** (à L² près) la forme
qui définit la projection, et elle ne vaut que pour l'élément **diagonal**.
Pour un hopping ⟨w₀|H|w_R⟩, les deux fonctions sont projetées sur des espaces
*différents* : plus de structure variationnelle, donc erreur du **premier**
ordre — attendue autour de ~12 %, pas ~1.5 %.

**Si ça se vérifie, la question §6-Q2 de la note 05 est tranchée** : le H¹ serait
trompeusement bon sur l'on-site et honnête sur les hoppings, ce qui justifierait
de réorienter la métrique. C'est **la** mesure à faire ensuite.

## 5. Provenance des données (source de confusion récurrente)

| Fichier | Date | nΦ | erreur H¹ | ζ max | Utilisable TB ? |
|---|---|---|---|---|---|
| `data/compressed_wannier_55_H1.json` | fév. 2025 | 10 | 5.59 % | 4.0 | ✅ |
| `data/compressed_wannier_55_H1_Ecut50.json` | fév. 2025 | 10 | 12.38 % | 1.3 | ✅ |
| `workflow/diag_outputs/diag_compressed_wannier.json` | 5 mai | 15 | 7.38 % | **5281.9** | ❌ §6 |

**Les deux fichiers de `data/` datent de février 2025 (ancien code) et ne sont
pas suivis par git.** Ce sont eux qui ont servi aux runs ci-dessus. Le seul
compressé produit par le code « nouveau » (Tikhonov) est celui de
`diag_outputs/`, inexploitable pour les intégrales.

## 6. Le piège ζ → ∞, mesuré sur un observable

Intégrales analytiques (code corrigé) :

| Fichier | ζ max | ‖w_gto‖ | T_gto | T/‖w‖² |
|---|---|---|---|---|
| `compressed_wannier_55_H1.json` | 4.0 | 0.9989 | 1.4782 | 1.4814 |
| `diag_compressed_wannier.json` | **5281.9** | **59.16** | **8.55e7** | **24435** |

Mécanisme (pathologie 1 de la note 02, appliquée au cinétique) :

- La fonction à ζ = 5282 est normée à 1 **analytiquement**, mais son contenu
  Fourier est concentré vers |q|² ~ 4ζ ≈ 21000, très au-delà de Ecut = 15 : sa
  représentation plane-wave est **quasi nulle**.
- Le greedy optimise dans le produit scalaire **discret**, où cette fonction
  « ne coûte rien » : il lui attribue un gros coefficient sans pénalité.
- L'intégrale analytique, elle, voit la vraie fonction : T ~ ζ, d'où 10⁷–10⁸.

Ce n'est pas une perte : le run de mai était un run de *diagnostic*
(`diagnose_conditioning.jl`, tol volontairement stricte), fait pour exhiber la
pathologie. Il fournit ici la preuve la plus nette de sa gravité (10⁷ sur un
observable, contre cond(S) = 2e7 dans la note 02).

## 7. État des remèdes au conditionnement (vérifié dans le code)

| Remède | Où | Statut | Ce qu'il traite |
|---|---|---|---|
| Barrière `ζ < ζ_min → Inf` | `compression.jl:14,67` | ✅ | Borne **inférieure** seulement (éviter les NaN) |
| Optimiseur du greedy | `compression.jl:17` | `ConjugateGradient()` **non contraint** | — |
| Borne supérieure ζ_max | greedy | ❌ **inexistante** | rien n'empêche l'escape |
| Tikhonov (solve intérieur) | `_make_inner_solver` | ✅ | conditionnement de S, **pas** l'escape |
| Box sur log ζ (sigmoïde) | `joint_optim/parametrization.jl` (`f6fd1bc`) | ✅ **phase B seulement** | borne des deux côtés, par reparamétrisation |
| Cohérence produits scalaires | — | ❌ **non fait** | la cause racine |

- **`Fminbox` n'apparaît nulle part dans l'historique.** La box de la phase B est
  un changement de variable sigmoïde, pas de l'optimisation contrainte.
- **Tikhonov n'a pas réglé l'escape** : le run de mai *avait* Tikhonov et a
  quand même produit ζ = 5282. Les deux remèdes agissent à des endroits
  différents — Tikhonov sur le *solve*, la borne sur la *sélection*.
- **La phase B est immunisée par construction** (ζ borné par la sigmoïde) : sa
  sortie sera directement utilisable en TB, contrairement au greedy.

**Question ouverte** : Tikhonov + ζ_max suffit-il, ou faut-il la cohérence des
produits scalaires ? Non tranché. Expérience simple : relancer le greedy avec
ζ_max ~ Ecut/4 et comparer, à nombre de fonctions égal, aux 5.59 % de février.

## 8. Rappel — qu'est-ce que la phase B ?

L'alternative au greedy. Le greedy ajoute les fonctions **une par une**, chacune
localement optimale, ce qui accumule de la redondance (pathologie 2, note 02) et
plafonne vers ~8.6 % à 11 fonctions.

La phase B optimise **tous les paramètres de toutes les fonctions
simultanément** — spreads ζ, centres/rayons π-bond, coefficients polynomiaux —
par L-BFGS avec gradients ForwardDiff, l'erreur H¹ comme loss et les
coefficients linéaires éliminés par **variable projection**. ζ est encodé par
sigmoïde, donc borné par construction.

État : Tasks 1-9 / 12 faites. Restent Task 10 (K-restarts + persistance),
11 (smoke), 12 (comparaison baseline phase A). **Décision en suspens** :
(C+A) run scaled-down ForwardDiff vs (B) passage à Enzyme — le gradient
ForwardDiff coûte ~42 s/itération à 113 paramètres, soit ~12 h pour la baseline.
Détails : `RESUME-HERE.md`, `03-design-phase-B.md`, `04-plan-phase-B.md`.

## 9. D'où vient le potentiel en tight binding ? (question ouverte)

Formellement H_{ij}(R) = ⟨w_i | −½∇² + V | w_j(·−R)⟩ où V est le potentiel du
**système cible**. Pour le TBG ce serait donc le V_KS du TBG — dont l'obtention
demanderait un SCF sur ~11 000 atomes, c'est-à-dire le calcul que le modèle TB
est censé remplacer.

Trois issues classiques :

1. **Slater-Koster empirique** : on ne calcule jamais V, on paramétrise
   directement t(R) (distance/orientation), ajusté sur l'expérience ou de la
   DFT petite cellule. C'est le gros de la littérature TBG.
2. **Superposition de potentiels monocouche** : V_TBG ≈ V₁(r) + V₂(r), chaque
   Vₗ étant le V_KS d'une monocouche translatée/tournée. Approximation
   non-auto-cohérente standard : on néglige la redistribution de charge due au
   twist, justifiée par le faible couplage van der Waals. **C'est la voie
   accessible et le socle implicite de presque tous les TB atomistiques du TBG.**
3. **Auto-cohérent complet** : seulement aux grands angles commensurables —
   sert de **benchmark pour valider (2)**.

**Position retenue pour ce projet : (2), validée par (3)** sur un grand angle
commensurable (étape 3 de l'échelle §4 note 05).

Subtilités à traiter le moment venu :
- **Champ cristallin / double comptage** : ⟨w_i^{(1)}|V₂|w_j^{(1)}⟩ n'est pas nul,
  c'est la correction exercée par l'autre couche. Physique, à garder — c'est ce
  que l'approximation à deux centres de Slater-Koster jette.
- **Coût** : superposer ne demande pas de SCF, mais représenter V sur la grille
  moiré reste lourd. Chaque Vₗ étant périodique à l'échelle de *sa* monocouche,
  on peut travailler dans la base d'ondes planes de la monocouche (petite) avec
  phase/rotation adéquates, sans construire la grille moiré — proprement aux
  angles commensurables.

**Pour L-2b (monocouche), la question ne se pose pas** : V_KS est sans ambiguïté
celui de `Graphene().scf()`. La difficulté n'arrive qu'au passage bicouche.

### À clarifier en session dédiée

Le but de la compression est de rendre calculables les ~10⁵–10⁶ intégrales de la
cellule moiré (chacune devenant une formule fermée au lieu d'une somme sur ~10⁶
ondes planes). Éviter le SCF est une **autre** question, traitée par la
superposition. Les deux sont nécessaires et orthogonales ; aucune méthode
n'échappe à la seconde approximation, et celle-ci est plus fine et plus
validable que le paramétrage empirique de Slater-Koster. **À discuter en détail
avant toute implémentation bicouche**, pour que l'objectif soit explicite.

## 10. Ce qu'il reste

**Chaîne de validation :**
1. **Hoppings t(R≠0)** — la mesure décisive (§4). Machinerie déjà presque là :
   côté référence, translater = multiplier par une phase e^{−iG·R} en Fourier ;
   côté gaussien, `translate` existe dans `CompressedWannier.jl`, et
   `integral(Wc₁, Wc₂; type=...)` accepte deux Wanniers distincts. Cinétique et
   overlap sont à portée immédiate, intralayer comme « pseudo-interlayer »
   (translation rigide de 3.35 Å en z).
2. **L-2b — terme V_KS on-site** : câbler un `scfres` cohérent
   cellule↔supercellule. Complète H = −½∇² + V_KS.
3. **T-band** : bandes monocouche vs DFT (cônes de Dirac, v_F).
4. Puis bilayer AB/AA → TBG grand angle → angle magique (échelle §4 note 05).

**Compression :**
- Relancer un greedy Ecut=15 avec ζ_max ~ Ecut/4 → un compressé « nouveau code »
  utilisable en TB ; comparer aux 5.59 % de février.
- Trancher §7 (Tikhonov + ζ_max vs cohérence des produits scalaires).
- Versionner les compressés de `data/`, ou au moins tracer leurs paramètres de
  génération (deux sessions y ont perdu du temps).

**Transverse :**
- Étendre `_check_L²_dot_precision` aux fonctions **décentrées** — c'est son
  angle mort qui a laissé passer le bug de convention pendant toute la phase A.
