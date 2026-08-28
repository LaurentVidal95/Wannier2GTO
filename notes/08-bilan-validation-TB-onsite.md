# Validation TB on-site — bilan

Bilan du premier maillon de la chaîne de validation tight-binding : comparaison
des intégrales TB calculées avec le vrai Wannier plane-wave et avec le Wannier
gaussien-compressé. Spec : `06-design-validation-TB-onsite.md`.
Plan : `07-plan-validation-TB-onsite.md`.

## TL;DR

- `laplacian_julia` / `kinetic_julia` sont **natifs** : la chaîne TB un corps
  est entièrement GaIn-free. Seul restait ce trou, il est bouché.
- Validation croisée à **trois routes indépendantes** (forme fermée, expansion
  `overlap_julia`, grille FFT) : accord 7e-7 à 1.4e-4 à Ecut=15, 1e-12 à Ecut=30.
- Deux runs on-site (Ecut=15 et 50) : l'erreur **normalisée** sur l'intégrale
  cinétique est de **0.89 %** à Ecut=50 pour une erreur H¹ de 12.4 %.
  **La métrique H¹ est très conservative pour l'élément on-site.**
- ⚠️ **Piège identifié** : le compressé recalculé en mai (`diag_outputs/`)
  contient une gaussienne à ζ ≈ 5282 qui fait exploser les intégrales
  analytiques d'un facteur 10⁸. Inutilisable pour le TB. Voir §4.

## 1. Ce qui a été livré

| Fichier | Rôle |
|---|---|
| `src/integrals/julia_integrals.jl` | `_axis_laplacian_terms`, `laplacian_julia`, `kinetic_julia` |
| `src/integrals/integrals.jl` | `symb_to_integral` : `:laplacian`/`:kinetic` → natif |
| `src/tight_binding/hamiltonian_scalar_product.jl` | `compare_onsite_kinetic` |
| `test/test_integrals_laplacian.jl` | 54 tests (forme fermée, grille, symétrie, H¹, intégration) |
| `workflow/validate_tb_onsite.jl` | script de comparaison, `W2G_ECUT` paramétrable |

Commits : `69104ca` → `eec1620`.

### Le principe

∇²X₂ reste une gaussienne-polynôme de **mêmes** ζ et R, donc tout se ramène à
des appels `overlap_julia` : pas de fonction de Boys, pas de McMurchie-Davidson,
et c'est **ForwardDiff-compatible**.

Convention de signe (imposée par `Hˢ_overlap(Ms; s=1)`, `Hs_scalar_prods.jl:32`) :

    laplacian_julia(X₁,X₂) = ⟨∇X₁,∇X₂⟩ = −⟨X₁,∇²X₂⟩     (forme positive)
    kinetic_julia          = ½ laplacian_julia

**Effet de bord utile** : `Hˢ_overlap(Ms; s=1)` (overlap H¹ analytique) tourne
désormais sans GaIn — potentiellement exploitable pour accélérer la phase B.

### Deux corrections par rapport au plan initial

1. **La forme fermée (2L+3)ζ du plan était fausse.** `dz2` (exposants (0,0,2))
   donne 1.3 au lieu de 2.1 à ζ=0.3. Redérivation par les moments 1D :

       ⟨∇g,∇g⟩ / ‖g‖² = ζ · Σᵢ (4nᵢ−1)/(2nᵢ−1)

   (2L+3)ζ n'est valable que si **tous** les exposants valent 0 ou 1 ; z²
   mélange L=2 et L=0, d'où 13ζ/3. Le code était correct, l'hypothèse ne
   l'était pas. Bonus : cette forme exerce la branche `a≥2` que le test s-only
   ne touchait jamais.

2. **Un couple de test du plan était dégénéré** : `xz`×`yz` s'annule par parité
   → l'assertion portait sur 0 ≈ 0, elle ne testait rien. Remplacé.

## 2. Provenance des données (important — source de confusion)

| Fichier | Date | nΦ | erreur H¹ | ζ max | Utilisable TB ? |
|---|---|---|---|---|---|
| `data/compressed_wannier_55_H1.json` | fév. 2025 | 10 | 5.59 % | 4.0 | ✅ |
| `data/compressed_wannier_55_H1_Ecut50.json` | fév. 2025 | 10 | 12.38 % | 1.3 | ✅ |
| `workflow/diag_outputs/diag_compressed_wannier.json` | 5 mai | 15 | 7.38 % | **5281.9** | ❌ voir §4 |

**Les deux fichiers de `data/` datent de février 2025 (ancien code) et ne sont
pas suivis par git.** Ce sont eux qui ont servi aux deux runs ci-dessous.
Le seul compressé produit par le code « nouveau » (Tikhonov) est celui de
`diag_outputs/`, et il est inexploitable pour les intégrales.

## 3. Résultats on-site

Setup commun : d = 10 Å, kgrid = [5,5,1], 10 fonctions de base.

| | Ecut=15 | Ecut=50 |
|---|---|---|
| erreur H¹ (compression) | 5.59 % | 12.38 % |
| ‖w_gto‖ | 0.9686 | 0.9793 |
| **T_ref** (référence grille) | **1.4721** | **1.5051** |
| T_gto (analytique) | 1.4261 | 1.4308 |
| erreur brute | 3.13 % | 4.94 % |
| **erreur normalisée** | **3.26 %** | **0.89 %** |

### Lecture

**La référence à Ecut=15 est sous-convergée.** T_ref passe de 1.4721 à 1.5051
(2.2 %) entre les deux runs. L'intégrale cinétique pèse les hautes fréquences
en |G|², donc une grille grossière tronque précisément ce qui contribue le
plus. **Toute conclusion tirée du seul run Ecut=15 est suspecte.**

En particulier, l'« inversion de signe » observée à Ecut=15 (le gaussien
normalisé sur-estimait : 1.5201 vs 1.4721), que j'avais d'abord lue comme la
signature physique d'une queue diffuse sous-représentée, **ne survit pas** à
Ecut=50 (1.4917 vs 1.5051, sous-estimation légère). C'était pour l'essentiel un
artefact de comparaison contre une référence tronquée.

Ce qui est solide, en revanche :

- **Erreur normalisée 0.89 % pour une erreur H¹ de 12.38 %** : facteur 14.
  Le H¹ est **très conservatif** pour l'élément on-site.
- **L'erreur brute est presque entièrement un déficit de norme** :
  ‖w_gto‖² = 0.959 et 1.4308/0.959 = 1.4917. L'erreur de *forme* est minuscule ;
  le compressé rate ~4 % de poids L², mais ce qu'il capture a la bonne courbure.

**Nuance** : l'on-site est le cas facile (intégrale locale, dominée par le
cœur). Ceci ne tranche pas la question §6-Q2 de la note 05 — ça dit que
**le problème n'est pas là**. Le hopping interlayer dépend de la queue p_z
hors-plan : c'est là qu'il faudra mesurer.

## 4. Le piège ζ → ∞, mesuré sur un observable

Intégrales analytiques calculées sur les deux compressés :

| Fichier | ‖w_gto‖ | T_gto | T/‖w‖² |
|---|---|---|---|
| `compressed_wannier_55_H1.json` (fév.) | 0.9686 | 1.4261 | 1.5201 |
| `diag_compressed_wannier.json` (mai) | **63.95** | **1.089e8** | **26629** |

Mécanisme (pathologie 1 de la note 02, appliquée au cinétique) :

- La fonction à ζ = 5282 est normée à 1 **analytiquement**, mais son contenu
  Fourier est concentré vers |q|² ~ 4ζ ≈ 21000, très au-delà de Ecut = 15 :
  sa représentation plane-wave est **quasi nulle**.
- Le greedy optimise dans le produit scalaire **discret**, où cette fonction
  « ne coûte rien » : il lui attribue un gros coefficient sans pénalité.
- L'intégrale cinétique analytique, elle, voit la vraie fonction : T ~ ζ,
  d'où 10⁸.

La référence-grille ne verrait rien de cette fonction, l'analytique la verrait
en entier : **la comparaison n'aurait aucun sens.**

Ce n'est pas une perte : le run de mai était un run de *diagnostic*
(`diagnose_conditioning.jl`, tol volontairement stricte « pour pousser
au-delà »), fait pour exhiber la pathologie. Il a rempli son rôle — et fournit
ici la preuve la plus nette de sa gravité (10⁸ sur un observable, contre
cond(S) = 2e7 dans la note 02).

## 5. État des remèdes — mise à plat (vérifié dans le code)

| Remède | Où | Statut | Ce qu'il traite |
|---|---|---|---|
| Barrière `ζ < ζ_min → Inf` | `compression.jl:14,67` | ✅ en place | Borne **inférieure** seulement (éviter les NaN) |
| Optimiseur du greedy | `compression.jl:17` | `ConjugateGradient()` **non contraint** | — |
| Borne supérieure ζ_max | greedy | ❌ **inexistante** | rien n'empêche l'escape |
| Tikhonov (solve intérieur) | `_make_inner_solver` | ✅ en place | conditionnement de S, **pas** l'escape |
| Box sur log ζ (sigmoïde) | `joint_optim/parametrization.jl` (`f6fd1bc`) | ✅ **phase B uniquement** | borne des deux côtés, par reparamétrisation |
| Cohérence produits scalaires | — | ❌ **non fait** | la cause racine |

Points à retenir :

- **`Fminbox` n'apparaît nulle part dans l'historique.** La box de la phase B
  est faite par changement de variable sigmoïde, pas par optimisation
  contrainte. (Le souvenir d'une « box optim avant Tikhonov » confond
  probablement les deux.)
- **Tikhonov n'a pas réglé l'escape** : le run de mai *avait* Tikhonov et a
  quand même produit ζ = 5282. Les deux remèdes agissent à des endroits
  différents — Tikhonov sur le *solve*, la borne sur la *sélection*.
- **La phase B est immunisée par construction** (ζ borné par la sigmoïde), donc
  sa sortie sera directement utilisable pour le TB, contrairement au greedy.
- Le correctif **structurel** reste à faire : normaliser dans le produit
  scalaire *discret* (ou tout passer en continu). ζ_max masque le mécanisme,
  la cohérence le supprime.

**Question ouverte** : Tikhonov + ζ_max suffit-il en pratique, ou faut-il la
cohérence des produits scalaires ? Non tranché. Une expérience simple :
relancer le greedy avec ζ_max ~ Ecut/4 et comparer, à nombre de fonctions égal,
avec les 5.59 % de février.

## 6. Rappel — qu'est-ce que la phase B ?

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

## 7. Suite

**Chaîne de validation :**
- **L-2b — terme V_KS** : seul point non trivial restant, câbler un `scfres`
  cohérent cellule↔supercellule (`Graphene().scf()` le recalcule). Complète
  H = −½∇² + V_KS on-site.
- **Hoppings t(R≠0)** : machinerie de translation (`translate` existe déjà dans
  `CompressedWannier.jl`) + shells de voisins. **C'est là que se joue la
  question de la métrique**, pas sur l'on-site.
- **T-band** : bandes monocouche vs DFT (cônes de Dirac, v_F).
- Puis bilayer AB/AA → TBG grand angle → angle magique (échelle §4 note 05).

**Compression :**
- Relancer un greedy Ecut=15 avec ζ_max ~ Ecut/4 pour obtenir enfin un
  compressé « nouveau code » utilisable en TB, et comparer aux 5.59 % de février.
- Trancher la question §5 (Tikhonov + ζ_max vs cohérence des produits scalaires).
- Envisager de versionner les compressés de `data/` (aujourd'hui non suivis),
  ou au moins de tracer leurs paramètres de génération.
