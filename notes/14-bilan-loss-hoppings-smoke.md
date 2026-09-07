# Bilan — smoke Ecut 15 de la loss hopping-ciblée (Jean Zay, 4 sept. 2026)

Exécution du plan [13](13-plan-loss-hoppings.md) Task 8, design
[12](12-design-loss-hoppings.md). Logs et script SLURM :
`workflow/joint_outputs/run_JeanZay_040926/` (5 jobs array, un par μ, 8 cœurs,
walltime 8 h, `--pkgimages=no`, N=10, K=3, 30 itérations, Ecut 15).

## 1. Résultats bruts

| μ (=ν) | état | H¹ rel. | C1 relerr T(a₁) | C2 signes val. | C3 \|S(a₁)\| | C4 H¹/H¹(μ=0) |
|---|---|---|---|---|---|---|
| 0 | fini (3 restarts) | **26.76 %** | FAIL (117 %) | FAIL (2 faux) | FAIL (1.58e-2) | 1.00× |
| 1 | fini (3 restarts) | 30.32 % | FAIL (185 %) | FAIL (3 faux) | **PASS (3.85e-4)** | 1.13× PASS |
| 10 | timeout 8 h | restart 1 fini (loss 0.154), restart 2 à l'itération 17 | — | — | — | — |
| 100 | timeout 8 h | restart 1 à l'itération 26 | — | — | — | — |
| 1000 | timeout 8 h | restart 1 fini (loss 0.438, 24 520 s), restart 2 à l'itération 1 | — | — | — | — |

Coût : ~105-160 s/itération à μ=0, ~155-240 s/itération à μ>0 (+50 %, le
terme de pénalité est loin d'être marginal — estimation du plan 13 fausse).
Un restart de 30 itérations ≈ 1h15 (μ=0) à 2-3 h (μ>0).

**Stalls** : à μ ≥ 1 uniquement, des itérations de 4 500-5 000 s (≈ 60-90
évaluations de line search) au terme desquelles la loss est *inchangée ou en
légère hausse* (+3e-7 relatif) — μ=1 : 1 stall ; μ=10 : 4 ; μ=100 : 5 ;
μ=1000 : 4. Jamais à μ=0. Ce sont ces stalls qui ont consommé le walltime.

## 2. Le bug : gradient de la pénalité faux sur les spreads

Vérification locale par différences finies (scratch, cas minuscule
N_centered=1, N_pibond=1, μ=ν=10) :

| terme | max écart AD vs FD | paramètres fautifs |
|---|---|---|
| H¹ | 1.5e-6 | aucun |
| pénalité hoppings (avant fix) | **1.0 (100 %)** | p1 et p9 = **les deux spreads ζ** |
| pénalité hoppings (après fix) | 5.0e-8 | aucun |

Cause racine : [GaussianPolynomials.jl:42](../src/basis_functions/GaussianPolynomials.jl)
calculait le préfacteur de normalisation via `filter_dual` — les partielles
étaient **délibérément supprimées** avant `analytic_norm` (contournement de
l'époque où cette norme passait par GaIn, incompatible avec les Duals). Le
chemin H¹ n'en souffrait pas (renormalisation sur la grille, en Dual, dans
`slow_fourier_transform_supercell`) ; mais le chemin analytique des hoppings
(`integral`) lit directement les coefficients normalisés → ∂/∂ζ faux. Fix :
suppression du filtre (les intégrales sont natives depuis longtemps) ; test de
régression FD ajouté dans `test/test_loss_hoppings.jl` ; 7 fichiers de tests
verts.

Conséquences : **tous les résultats μ>0 de ce smoke sont invalides** (L-BFGS
combattait un gradient incohérent — d'où les stalls HZ). Le PASS du critère 3
à μ=1 est directionnellement encourageant (la pénalité ν pousse bien S(a₁)
vers 0 : 1.58e-2 → 3.85e-4) mais ne sera acté qu'après relance.

Leçon de méthode : le plan 13 testait le gradient "fini et non nul", pas
"juste" — un FD est obligatoire pour tout nouveau terme de loss, même quand le
backend AD est validé sur un autre terme.

## 3. Le problème de fond : la base de comparaison

Même le témoin μ=0 (loss H¹ pure, gradient correct) est à **26.8 % après 30
itérations, loin des 5.6 % du greedy** — loss encore en pleine descente, T(a₁)
de signe faux à μ=0. La question du design (la pénalité répare-t-elle les
hoppings à forme égale ?) n'est pas encore *posable* : on compare des
pénalités sur une fonction qui n'est pas ajustée.

C'est l'écart assumé du plan 13 (init `init_params` + K restarts au lieu du
greedy) qui se réalise. Et l'encodage greedy→`JointLayout` n'est **pas
possible tel quel** (vérifié sur `data/compressed_wannier_55_H1.json`) : chaque
fonction greedy est une *contraction* de 5-7 SAGTOs de **spreads différents**
(ex. Φ1 : spreads 0.66/1.21/1.17/0.29/0.36…), alors qu'une fonction
`JointLayout` a **un seul ζ** partagé par ses groupes polynomiaux. Représenter
le greedy dans le layout demanderait ~30-40 fonctions (≈400 paramètres).

## 4. Le problème de coût

~0.5 s par évaluation de loss × 113 partielles ForwardDiff ≈ 55-80 s par
gradient, 2-3 évaluations par itération. À ce tarif, les runs sérieux du
design (Ecut 50, N=15, 200 itérations, K=5) sont hors de portée, même en
parallèle par μ.

Levier structurel identifié : **le jacobien est bloc-creux** — chaque
paramètre n'affecte que *sa* fonction. En mettant en cache les évaluations
Fourier des fonctions inactives (Float64) et en ne passant en Dual que le bloc
actif (≈10 partielles), le gradient coûte ≈ N_blocs × (1/N_blocs) × 10
évaluations ≈ 10 évaluations de loss au lieu de 113 → **~10×**, sans nouvelle
math ni changement de backend. Même logique pour `gto_hoppings` (seules les
paires impliquant le bloc actif à recalculer).

## 5. Décisions et suite

1. **Fix gradient : fait**, à committer avec le test de régression et
   `W2G_K` (nombre de restarts par env — K=1 suffit pour un scan apparié, les
   seeds étant partagées).
2. **Relance cheap de validation** (cluster) : μ ∈ {0, 1, 10}, K=1, 30
   itérations, 3 jobs parallèles (~1h15-3 h chacun). Objectif : confirmer la
   disparition des stalls sur le vrai problème et lire la tendance ν/S(a₁) —
   pas encore le verdict des critères (cf. §3).
3. **Design 15 — gradient bloc-cache** : prérequis de tout run sérieux.
   Petit design, gain attendu ~10×.
4. **Design 16 — init greedy** : nouveau layout "à spreads par groupe"
   (structure du greedy : 7 groupes centrés, 5 groupes π-bond × 3 copies D3)
   pour partir de la solution greedy à 5.6 % et ne raffiner que la queue. Sans
   ça, la question du design 12 ne peut pas être tranchée.
5. Ordre : 2 (lance et oublie) → 3 → 4 → vrai smoke → runs Ecut 50.

Le smoke a coûté ~40 h·job de cluster et n'a pas tranché le design — mais il a
débusqué un bug de gradient latent qui aurait pollué *toute* la suite, et
mesuré les deux verrous (init, coût) qu'il fallait lever de toute façon.
