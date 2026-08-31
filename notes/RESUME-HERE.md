# Reprise — état au 31 août 2026 (soir)

## Où on en est

**Loss hopping-ciblée : implémentée et committée en entier** (design
[12](12-design-loss-hoppings.md) + biblio [12b](12b-biblio-loss-hoppings.md),
plan [13](13-plan-loss-hoppings.md), HEAD `fc7018c`). 131 tests verts.

- `reference_hopping` / `HoppingTargets` / précalcul des cibles (seed 1234,
  `data/` gitignoré → régénérer via
  `W2G_ECUT=15 julia --project=. workflow/precompute_hopping_targets.jl`) ;
- `gto_hoppings` + `hopping_penalty` branchés dans `joint_loss`
  (ForwardDiff vert sans relâchement de types) ;
- `run_joint_optim` (ex-Task 10 phase B) avec passthrough μ/ν ;
- évaluation par critère (PASS/FAIL indépendants, demande de Laurent) ;
- `workflow/joint_compression_hoppings.jl` : scan μ apparié (μ ∈ {0,1,10,100,1000},
  ν=μ, MASTER_SEED partagé = paires contrôlées — ne pas "corriger" ça).

## En attente : les runs (côté Laurent, sur cluster)

Le smoke local a été stoppé (trop long laptop). Laurent lance sur cluster :
précalcul cibles puis smoke Ecut 15 (`W2G_MAX_ITER=30`), puis si concluant les
runs sérieux Ecut 50 (`W2G_MAX_ITER=200`, monter N_CENTERED/N_PIBOND/K dans le
script). Sanity précalcul : intra a1 T_ref ≈ 1.30e-2 (Ecut 15) / 3.89e-2
(Ecut 50) ; inter AA S_ref ≈ −4.95e-2 / −5.52e-2.

## Pour reprendre

1. Demander les logs (`workflow/diag_outputs/joint_hoppings_*.log`).
2. Écrire `notes/14-bilan-loss-hoppings-smoke.md` : tableau μ ↦ (H¹, critères
   1-4), quel critère mord, trade-off, go/no-go Ecut 50 (plan 13, Task 8).
3. Critères (design 12 §5) : (1) relerr T(a₁) ≤ 5 % ; (2) signes corrects
   validation si |T_ref| > 5e-4 ; (3) |S(a₁)| ≤ 1e-3 ; (4) H¹ ≤ 1.5× jumeau
   μ=0. Attention : si μ=0 n'atteint pas ~5.6 % H¹ (niveau greedy), l'écart
   assumé du plan 13 (init `init_params` au lieu de greedy) devient le suspect
   n°1 → mini-design encodage greedy→JointLayout.
4. Ensuite : approche B (t(R) complet, scfres L-2b — piste
  `../../TBG/TwistedBilayerGraphene.jl`, chapitre de thèse dispo sur demande)
  en validation finale ancrée littérature (t₁ ≈ −2.7 eV, t₂/t₁).

## Repères techniques

- Tests : `test_hoppings.jl` (18), `test_hopping_targets.jl` (30),
  `test_loss_hoppings.jl` (11), `test_loss.jl` (10),
  `test_integrals_laplacian.jl` (62) — tous fichier par fichier.
- v2 documentées non implémentées (design 12 §2) : lagrangien augmenté ;
  formulation riemannienne "shift-orthonormalité" — décision Laurent (sa
  thèse partie 1 : Stiefel/Flag).
- Contexte scientifique : bilan [11](11-bilan-validation-hoppings.md) (T ordre
  1, signes faux, orthonormalité violée) ; position actée : T doit être bon
  *en soi*, pas de compensation V_KS (non transportable au TBG).
