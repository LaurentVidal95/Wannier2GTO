# Reprise — état au 5 mai 2026

## Où on en est

Phase B en cours d'exécution via `superpowers:subagent-driven-development`.

**9 tâches sur 12 faites** :
- ✅ Task 1: deps & module skeleton (`a8218bf`)
- ✅ Task 2: `JointLayout` + flat ↔ structured (`c0795097`)
- ✅ Task 3: sigmoid encoding pour log ζ (`f6fd1bc`)
- ✅ Task 4: `params_to_basis_functions` avec D3 (`f3bc173`)
- ✅ Task 5: `init_params` heuristique (`0a0b284`)
- ✅ Task 6: `joint_inner_solve` Tikhonov (`bc44c11`)
- ✅ Task 7: `joint_loss` H¹ via variable projection (`833cd94`)
- ✅ Task 8 v2: gradient ForwardDiff validé (`c9b5d19`)
- ✅ Task 9: `run_lbfgs_once` (`4087c4f`)
- ⏳ Task 10: K-restarts + persistance — **à faire**
- ⏳ Task 11: smoke workflow — **à faire**
- ⏳ Task 12: comparaison baseline phase A — **à faire**

## Décision en suspens

Le benchmark Task 9 a montré le coût du gradient ForwardDiff (~7.5 s / iter à
20 paramètres, extrapolé à ~42 s / iter à 113 paramètres) :

| Run | Coût estimé |
|---|---|
| Smoke (N=10, K=3, 50 iter, Ecut=15) | ~30 min |
| Baseline (N=15, K=5, 200 iter, Ecut=15) | ~12 h |
| Ecut=30 baseline | ~60 h (cluster) |

J'avais proposé trois options :

- **(A) On continue ForwardDiff** : Tasks 10-11, smoke local, baseline overnight ou cluster
- **(B) On tente Enzyme maintenant** : potentiellement 5× plus rapide mais risque (cf. échec Zygote)
- **(C) Run scaled-down d'abord** : N=10, K=3, max_iter=100 → ~2.5 h laptop pour un proof-of-concept

Et recommandé **(C) puis (A)**. **Validation utilisateur en attente** au moment de la pause.

## Pour reprendre

1. Reouvre la session : `claude --continue` depuis le répertoire du projet
2. Ma première phrase à la reprise sera de te re-poser la question (C+A) ou (B) — tu m'as dit "Mince, je vais devoir m'arrêter là momentanément" sans avoir répondu
3. Si tu valides (C+A) : je dispatche Task 10 (orchestration K-restarts), puis Task 11 (smoke), puis on lance le smoke et on analyse

## Repères techniques

- HEAD courant : `4087c4f Phase B: single-run L-BFGS with ForwardDiff gradients`
- Tests : `julia --project=. test/test_loss.jl` (10 tests verts), `julia --project=. test/test_init.jl` (3 testsets verts), `julia --project=. test/test_parametrization.jl` (5 testsets verts)
- Smoke pipeline : `julia --project=. -e '...'` + `include("test/test_loss.jl")` charge le `_BASIS_SC` (DFTK supercell, ~30 s) puis l'optim peut tourner
- Limites identifiées :
  - Cold start ForwardDiff ≈ 10 min de compilation
  - Warm gradient ≈ 7.5 s à 20 params, ~42 s à 113 params
  - 30% GC time, 207 GB allocs / 10-iter run → ForwardDiff alloue beaucoup

## Spécifications de référence

- `notes/03-design-phase-B.md` — design (architecture, hyperparams, success criteria)
- `notes/04-plan-phase-B.md` — plan d'implémentation 12 tâches
- `notes/02-bilan-phase-A.md` — résultats phase A (greedy 8.6% à 11 fonctions)
