# Phase A — Diagnostic et remédiation du conditionnement de $S$

## 1. Symptôme

D'après thèse 5.4.2.3 :
> *Adding more basis functions to $\Phi^{SA}$ resulted in the conditioning of
> the overlap matrix $S$ to blow-up. This main limitation of our proof-of-concept
> code, which hindered the compression procedure, can be remedied for example
> by minimizing the condition number of $S$ in the compression procedure, or
> by applying other standard methods to cure ill-conditioned basis sets.*

État connu : avec $n_x = n_z = 5$, $\Phi^{\mathrm{SA}} \in \mathcal{B}_{5,5}$ de
10 fonctions, 155 paramètres, $\|w_z - \Pi w_z\|_{H^1} \simeq 0.12$. Au-delà,
$S$ explose.

## 2. Causes probables (à vérifier)

Hypothèses physiques sur l'origine du mauvais conditionnement :

- **(H1) Quasi-degénérescence en spread** : deux SAGTOs centrées au même
  endroit avec des $\zeta$ proches deviennent quasi-colinéaires.
  C'est la cause classique en chimie quantique avec des bases $sp$ étendues
  (cf. Lehtola 2019b).
- **(H2) Quasi-redondance entre symétries D3 vs non-D3** : l'alternance
  D3_sym / non-D3_sym dans `compression.jl:128–129` ajoute des SAGTOs
  centrées sur les liaisons $\pi$ ; combinées à 3 rotations $D_3$, elles
  peuvent recouvrir partiellement les SAGTOs centrées.
- **(H3) Polynômes proches en orthogonalité numérique** : à grand ordre,
  les polynômes de la table 5.2 deviennent oscillants et leurs intégrales
  numériques sont sujettes à perte de précision.

À vérifier par instrumentation (étape 3 ci-dessous).

## 3. Instrumentation à ajouter (avant tout fix)

But : poser un diagnostic chiffré, pas patcher à l'aveugle.

### 3.1. Logger à chaque itération du greedy

À ajouter dans `compress_graphene_pz_wannier` :

```julia
# After computing S = Hˢ_overlap(...)
σ = svdvals(S)
@info "iter=$n_iter cond(S)=$(σ[1]/σ[end]) σ_min=$(σ[end]) σ_max=$(σ[1]) n_AOs=$(length(SAGTOs))"
```

Idéalement, stocker dans `info` (struct callback) :
- `cond(S)`, $\sigma_{\min}(S)$, $\sigma_{\max}(S)$
- distribution des spreads $\zeta_i$ et des centres $\alpha_i$
- norme du résidu et erreur courante

### 3.2. Identifier *quelle* SAGTO ajoutée fait basculer le conditionnement

Tracer cond$(S_n)$ en fonction de $n$ (croissance attendue exponentielle quand
on entre en quasi-degénérescence). Identifier les itérations critiques.

### 3.3. Vérification croisée (H1)

Pour chaque paire $(i, j)$ à l'itération critique : afficher
$|\langle \phi_i, \phi_j \rangle| / \sqrt{\|\phi_i\| \|\phi_j\|}$.
Quasi-1 ⇒ paire colinéaire identifiée.

## 4. Remédiations à tester (par ordre de complexité)

### Option A1 — Régularisation de Tikhonov dans la résolution intérieure

```julia
# Replace S \ Γ by:
ε = 1e-8 * tr(S) / size(S, 1)  # scale-invariant
optimal_coeffs = (S + ε*I) \ Γ
```

Avantage : 1 ligne de diff, immédiat.
Inconvénient : biaise la projection ; à utiliser comme *baseline*, pas comme solution finale.

### Option A2 — Filtrage par décomposition propre de $S$ avec seuil

```julia
F = eigen(Hermitian(S))
keep = F.values .> tol_sval * maximum(F.values)
# Project onto kept subspace, solve, re-expand
```

Avantage : projette dans un sous-espace numériquement stable, sans biais de Tikhonov.
Inconvénient : la SAGTO ajoutée à l'itération courante peut être filtrée → potentielle perte d'expressivité ; signal qu'on a ajouté une fonction redondante.

### Option A3 — Cholesky pivotée (Lehtola 2019b)

Approche standard pour bases gaussiennes surcomplétées. Une factorisation
$P^T S P = L L^T$ avec pivotage et seuil tronque proprement les directions
quasi-nulles.

Avantage : méthode standard, théorie propre, robuste.
Inconvénient : un peu plus de code (orthogonalisation des SAGTOs *avant* le
calcul des coefficients de combinaison).

Implémentation Julia : `LinearAlgebra.cholesky(Hermitian(S), RowMaximum(), check=false)` (Julia 1.9+).

### Option A4 — Pénaliser cond$(S)$ dans l'objectif extérieur

Modifier $f$ dans `compression.jl:58–71` :
```julia
err = ... # current Hˢ error
penalty = log(cond(S))  # or smooth surrogate
return err + λ * penalty
```

Avantage : fait par construction le tradeoff précision/conditionnement, suggéré dans la thèse.
Inconvénient : ajoute un hyperparamètre $\lambda$ ; cond$(S)$ n'est pas trivialement
différentiable en mode auto (utiliser `logabsdet` + `tr` plutôt).

## 5. Ordre d'attaque proposé

1. **Instrumentation (§3.1, §3.2)** — 30 min
2. **Reproduire le bug** : faire tourner la compression existante jusqu'à blow-up,
   logger pour valider H1/H2/H3 — 1h
3. **A1 (Tikhonov)** comme sanity check : confirme que c'est bien
   le conditionnement qui bloque et pas un autre bug — 30 min
4. **A3 (Cholesky pivotée)** comme remédiation propre — 2–3h
5. **A4 (régul cond$(S)$)** si nécessaire pour aller plus loin — différé

## 6. Critère de succès de la phase A

- Pouvoir faire tourner la compression jusqu'à au moins **20 fonctions de base**
  sans explosion numérique
- Erreur $H^1$ relative $\le 5\%$ (à comparer aux 12% actuels avec 10 fonctions)
- cond$(S)$ resté $\le 10^{10}$ tout au long du greedy

## 7. Hors scope de A

- Migration vers Flux.jl
- AD reverse-mode (Zygote/Enzyme) — peut être tenté en bonus si simple
- Reformulation conjointe (= phase B)
- Calculs aval ERI / tight-binding TBG

## 8. Livrables

- Patch sur `src/compression/compression.jl` + `src/common/Hs_scalar_prods.jl`
  si overlap calculé là
- Un script `workflow/diagnose_conditioning.jl` qui reproduit et visualise
  cond$(S)$ vs itération
- Mise à jour de cette note avec les conclusions chiffrées
