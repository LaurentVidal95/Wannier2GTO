# Design — validation on-site de la chaîne TB un corps

Date : 2026-06-24. Spec issu du brainstorm de reprise (thread « chaîne de
validation tight-binding », cf. `05-TB-landscape-et-validation.md`).

## 0. Objectif et scope

Clore (au premier maillon) la chaîne de validation TB un corps en **comparant
les intégrales TB calculées avec le vrai Wannier plane-wave et celles calculées
avec le Wannier gaussien-compressé**.

- **Hamiltonien** : un corps `H = −½∇² + V_KS` (décision note 05 §2).
- **Échelle de cette session** : on-site seulement, `⟨w|H|w⟩` (monocouche,
  Wannier p_z unique). Les hoppings `t(R)` (R≠0) et la reconstruction de bandes
  sont **à terme**, hors scope ici (échelle §4 de la note 05, étapes 2+).
- **Métrique** : on **garde** le Wannier compressé H¹ existant (phase A) tel
  quel et on **mesure** l'erreur relative sur les intégrales physiques. La
  comparaison est l'expérience qui informera plus tard un éventuel pivot de la
  loss de phase B vers une métrique « hopping » — pas de redéfinition de loss
  maintenant (décision brainstorm Q2).

## 1. État du code pertinent (rappel)

- `src/integrals/julia_integrals.jl` — `overlap_julia` natif (exact,
  ForwardDiff-compatible). Convention : `X(r) = x^{nx} y^{ny} z^{nz}
  exp(−ζ‖r−R‖²)`, **monômes en coordonnées absolues**, gaussienne centrée en R.
- `src/integrals/integrals.jl` — dispatch `symb_to_integral`. `:overlap` →
  natif ; `:laplacian`, `:kinetic` → `GaIn` (mort sur cette machine).
- `src/common/Hs_scalar_prods.jl` —
  - `Hˢ_dot(basis, ψ1, ψ2; s)` : poids de Sobolev `(1+|G|²)^s` en Fourier.
  - `Hˢ_overlap(Ms; s=1)` (ligne 32) fait `overlap + integral(:laplacian)` →
    **fixe la convention de signe** (voir §2).
  - `_check_L²_dot_precision` : patron de l'auto-validation analytique vs grille.
- `src/tight_binding/hamiltonian_scalar_product.jl` — `hamiltonian_scalar_prod`
  somme `[:kinetic, :PBE_local_potential]` ; `potential_scalar_prod` gère
  `GaussianPolynomial` / `BasisFunction` / `CompressedWannier`.
- Un `CompressedWannier` porte **les deux** représentations sur le **même**
  `basis_supercell` : `Wc.wannier` (vrai Wannier en Fourier) et
  `Wc.basis_functions` / `Wc.coefficients` (gaussien).

## 2. L-1 · `laplacian_julia` — maths, convention, signature

### Identité

Pour `X₂ = x^a y^b z^c · g`, `g = exp(−ζ‖r−R‖²)`, le Laplacien reste une
gaussienne-polynôme de **mêmes** ζ et R. Par direction (x ; idem y, z) :

```
∂²ₓ X₂ = [ a(a−1) x^{a−2}
           − 2ζ x^a
           − 4ζ a (x−Rₓ) x^{a−1}
           + 4ζ² (x−Rₓ)² x^a ] · g
```

`∇²X₂ = ∂²ₓX₂ + ∂²ᵧX₂ + ∂²_zX₂`. Chaque `(x−Rₓ)^k x^m` se réexpand en monômes en
coordonnées absolues (binôme), donc `∇²X₂` est une **somme finie de SAGTOs** de
spread ζ et centre R. D'où

```
⟨X₁, ∇²X₂⟩ = Σₖ coeffₖ · overlap_julia(X₁, monômeₖ)
```

— uniquement des appels `overlap_julia` : exact, sans Boys ni
McMurchie-Davidson, **ForwardDiff-compatible**.

### Convention de signe (NON négociable)

`Hˢ_overlap(Ms; s=1)` ajoute `integral(:laplacian)` à l'overlap L² pour obtenir
l'overlap H¹. Le produit scalaire H¹ étant `⟨X₁,X₂⟩ + ⟨∇X₁,∇X₂⟩`, on a

```
laplacian_julia(X₁, X₂) := ⟨∇X₁, ∇X₂⟩ = −⟨X₁, ∇²X₂⟩   (forme positive)
kinetic_julia(X₁, X₂)   := ½ laplacian_julia(X₁, X₂)    = ⟨X₁, −½∇² X₂⟩
```

Implémentation : on calcule `⟨X₁,∇²X₂⟩` via l'identité ci-dessus, puis on
**renvoie son opposé**.

### Signature et câblage

Même signature scalaire que `overlap_julia` (drop-in pour `symb_to_integral`) :

```julia
laplacian_julia(ζ1, R1, nx1, ny1, nz1, ζ2, R2, nx2, ny2, nz2) -> Real
```

Brancher dans `src/integrals/integrals.jl` :
`symb_to_integral[:laplacian] = laplacian_julia`. Le `:kinetic` reste dérivé
(½ laplacian) ; à décider à l'implémentation : entrée dédiée `kinetic_julia` ou
facteur ½ au point d'usage. Les niveaux `BasisFunction` / `CompressedWannier`
héritent automatiquement via les boucles `integral(...)` existantes.

### Auto-validation contre la grille (sans GaIn)

`∇²` = multiplication par `−|G|²` en Fourier. Dans l'esprit de
`_check_L²_dot_precision` :

```
laplacian_julia(X₁,X₂)  ≈  dot( X₁(basis_sc), |G|² .* X₂(basis_sc) )   # ⟨∇X₁,∇X₂⟩
```

avec `|G|²` issu de `G_vectors_cart(basis_sc, only(basis_sc.kpoints))`. Tester à
ζ modéré (≲ Ecut/4) pour éviter la pathologie-1 d'escape (note 02).

## 3. L-2 · Harnais de comparaison on-site

Sur un `CompressedWannier` chargé, pour `H = −½∇² + V_KS` :

| Terme | Référence (vrai Wannier, grille) | Gaussien (compressé, analytique) |
|---|---|---|
| Cinétique `⟨w|−½∇²|w⟩` | `½ · real(dot(w_four, |G|² .* w_four))` | `½ · laplacian_julia` via coeffs (= `integral(Wc,Wc; type=:laplacian)/2`) |
| Potentiel `⟨w|V_KS|w⟩` | quadrature réelle `∫ w(r)² V_KS(r) dr` | `potential_scalar_prod(basis, V_KS_four, Wc, Wc)` (existe ✓) |

`w_four = Wc.wannier`. Sortie du harnais : pour chaque terme disponible,
`(valeur_réf, valeur_gauss, erreur_relative)`, plus rappel de `Wc.error_norm`
(H¹). Le gaussien on-site = `integral(Wc, Wc; type=:laplacian)` désormais natif.

### Séquencement (dé-risque)

1. **L-2a — terme cinétique** : autosuffisant, **zéro `scfres`**. Livré et
   validé avec L-1.
2. **L-2b — terme V_KS** : nécessite un `scfres` cohérent avec le
   `basis_supercell` du Wannier (subtilité cellule↔supercellule). Branché dans
   un second temps. `scfres` recalculable via `workflow/monolayer_graphene.jl`
   (`Graphene().scf()`, PBE).

## 4. Données

- Wannier compressé (gaussien) : `data/compressed_wannier_55_H1_Ecut50.json`.
- Vrai Wannier (référence) : `workflow/wannier_functions/wannier_pz_Ecut-50.json`.
- Ecut=50 : référence-grille la plus précise → test le plus propre de
  l'analytique (`laplacian_julia` est exact ; l'écart résiduel = résolution de
  grille). La comparaison à Ecut=50 sera **lancée sur cluster** par
  l'utilisateur.

## 5. Tests

`test/test_integrals_laplacian.jl` :

- **(a) analytique vs grille** : `laplacian_julia` ≈ `dot(X₁(basis), |G|² .*
  X₂(basis))` sur quelques paires de SAGTOs ciblées (ζ modéré).
- **(b) symétrie** : `laplacian_julia(X₁,X₂) == laplacian_julia(X₂,X₁)`.
- **(c) cas fermé** : gaussienne s pure normalisée, `⟨g,−∇²g⟩` a une forme
  analytique connue → vérifie l'amplitude absolue (pas seulement la cohérence
  interne grille).
- **(d) H¹ overlap** : `Hˢ_overlap(Ms; s=1)` tourne désormais sans GaIn et
  reste SPD sur un petit jeu de SAGTOs.

Exécution des tests : par l'utilisateur, en local
(`julia --project=. test/test_integrals_laplacian.jl`). Pas d'exécution ni
d'opération git autonome.

## 6. Hors scope (rappels pour plus tard)

- Hoppings `t(R)`, R≠0 : machinerie de translation (phase en Fourier côté réf,
  décalage de centres côté gaussien) + énumération des shells de voisins.
- Reconstruction de bandes (T-band), montée bilayer / TBG, termes deux corps.
- Pivot éventuel de la métrique de compression vers les hoppings (alimenté par
  les mesures de ce harnais).
