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
- Smoke Ecut=15 (10 fonctions) : erreur H¹ **5.59 %**, erreur sur l'intégrale
  cinétique on-site **3.13 %**.
- Signal intéressant : le compressé perd ~6 % de poids L² et ce qui reste est
  **trop courbé par unité de norme** — cohérent avec une queue diffuse
  sous-représentée (cf. drapeau §5 de la note 05).

## Ce qui a été livré

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

## Deux corrections par rapport au plan initial

1. **La forme fermée (2L+3)ζ du plan était fausse.** `dz2` (exposants (0,0,2))
   donne 1.3 au lieu de 2.1 à ζ=0.3. Redérivation par les moments 1D :

       ⟨∇g,∇g⟩ / ‖g‖² = ζ · Σᵢ (4nᵢ−1)/(2nᵢ−1)

   (2L+3)ζ n'est valable que si **tous** les exposants valent 0 ou 1 ; z²
   mélange L=2 et L=0, d'où 13ζ/3. Le code était correct, l'hypothèse ne
   l'était pas. Bonus : cette forme exerce la branche `a≥2` que le test s-only
   ne touchait jamais.

2. **Un couple de test du plan était dégénéré** : `xz`×`yz` s'annule par parité
   (intégrale en y impaire autour du centre commun) → l'assertion portait sur
   0 ≈ 0, elle ne testait rien. Remplacé par des paires non dégénérées.

## Résultat du smoke (Ecut=15, d=10 Å, kgrid=[5,5,1], 10 fonctions)

    compression error (H¹)    = 5.59 %
    ‖w‖ référence / gaussien  = 1.0000 / 0.9686

    brut :          T_ref = 1.4721    T_gto = 1.4261    → 3.13 %
    normalisé :     T_ref = 1.4721    T_gto = 1.5201    → 3.26 %

### Lecture (signaux, pas conclusions — un seul point, on-site, Ecut=15)

- L'erreur sur l'intégrale cinétique (3.1 %) est **plus petite** que l'erreur H¹
  (5.6 %) → la métrique H¹ est plutôt **conservative** ici. Premier élément de
  réponse à la question §6-Q2 de la note 05.
- **Le signe s'inverse à la normalisation** : brut le gaussien sous-estime
  (1.426 < 1.472), normalisé il sur-estime (1.520 > 1.472). Le compressé perd
  ~6 % de poids L² et ce qui reste est trop courbé par unité de norme —
  signature d'un greedy qui capture le cœur piqué et sous-représente la queue
  diffuse. Or **c'est la queue p_z hors-plan qui gouverne le hopping
  interlayer** (drapeau §5 note 05). Appui quantitatif préliminaire à
  l'hypothèse « la métrique pertinente n'est peut-être pas H¹ ».

## Prochaine action immédiate

Run cluster à Ecut=50 :

    W2G_ECUT=50 julia --project=. workflow/validate_tb_onsite.jl

(Le script asserte la cohérence Wannier/base — nombre de vecteurs G,
normalisation, finitude — avant tout calcul. Données requises :
`data/compressed_wannier_55_H1_Ecut50.json` et
`workflow/wannier_functions/wannier_pz_Ecut-50.json`, présentes.)

Point d'attention : vérifier que le compressé Ecut=50 a bien été produit avec
d=10 Å et kgrid=[5,5,1] ; sinon l'assert `length(w) == n_G` le signalera.

## Suite de la chaîne

- **L-2b — terme V_KS** : seul point non trivial restant, câbler un `scfres`
  cohérent cellule↔supercellule (`Graphene().scf()` le recalcule). Complète
  H = −½∇² + V_KS on-site.
- **Hoppings t(R≠0)** : machinerie de translation (phase en Fourier côté
  référence, décalage de centres côté gaussien — `translate` existe déjà dans
  `CompressedWannier.jl`) + énumération des shells de voisins.
- **T-band** : reconstruction de la structure de bandes monocouche, comparaison
  DFT (cônes de Dirac, v_F).
- Puis montée bilayer AB/AA → TBG grand angle → angle magique (échelle §4
  note 05).

## Thread parallèle : phase B

Intacte, indépendante. Tasks 10-12 restantes et décision en suspens
(C+A ForwardDiff vs B Enzyme) : voir `RESUME-HERE.md`.

Note : `Hˢ_overlap(Ms; s=1)` étant désormais natif et AD-compatible, il y a
peut-être là un gain à récupérer sur le coût du gradient (à évaluer).
