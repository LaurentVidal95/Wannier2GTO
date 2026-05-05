# Phase B — design : fit conjoint des SAGTOs

Note de spécification produite à l'issue du brainstorming phase B (mai 2026).
Objectif : remplacer la procédure greedy actuelle par un fit conjoint de tous
les paramètres pour dépasser le plafond identifié en phase A
(8.6% H¹ à 11 fonctions). Toutes les conclusions de phase A sont supposées
intégrées (voir `02-bilan-phase-A.md`).

## 1. Critère de succès

**Cible primaire** : matcher Bakhta+18, soit **≤ 5% H¹ relatif** sur la
projection de $w_z$, à nombre de fonctions comparable (~10–20). C'est ce qui
justifie l'investissement par rapport au greedy.

**Cible secondaire** (hors phase B mais oriente le design) : la précision
finale doit être suffisante pour que les calculs tight-binding aval (TBG)
soient utilisables. Si après la phase B les erreurs TB sont trop grandes,
on poussera vers une cible plus serrée (1–2% H¹) en phase ultérieure.

## 2. Architecture

Modèle paramétré :
$$w_z(\mathbf{r}) \approx \sum_{i=1}^{N} c_i\, \Phi_i(\mathbf{r}; \theta_i)$$
où chaque $\Phi_i$ est une SAGTO $A_2''$-symétrique (architecture A
validée pendant le brainstorming, équivalente à la structure du code
existant).

**Deux types de fonctions de base** :

| Type | Centre | Params libres | Total |
|------|--------|---------------|-------|
| Type 1 — *centré* | $\alpha = 0$ (forcé) | $\log\zeta$, 6 coefs poly $\{\lambda_{ij}\}$, $c_i$ | 8 |
| Type 2 — *π-bond* | $\alpha = r\,\hat{\mathbf{u}}_\pi$ avec $\alpha_z = 0$ forcé | $r$, $\log\zeta$, 4 coefs poly, $c_i$ | 7 |

Pour chaque fonction non-centrée, la symétrie $D_3$ est imposée par
construction via `enforce_D3_symmetry` : une seule paire $(r, \log\zeta)$
décrit trois groupes de SAGTOs aux centres rotationnellement équivalents.

Pour $N_\text{centered} + N_\text{pibond} = N = 15$ avec split 8/7 :
$$\dim\theta = 8\cdot 8 + 7\cdot 7 = 113 \text{ paramètres}$$

## 3. Frame d'optimisation

**Variable projection** : les coefficients linéaires $c_i$ sont éliminés à
chaque évaluation de la loss par closed-form solve Tikhonov-régularisé,
$$\mathbf{c}^\star(\theta) = (S(\theta) + \varepsilon I)^{-1}\,\Gamma(\theta),$$
où $S_{ij} = \langle\Phi_i, \Phi_j\rangle_{H^1}$,
$\Gamma_i = \langle w_z, \Phi_i\rangle_{H^1}$,
$\varepsilon = 10^{-8}$ (paramètre repris de phase A).

L'optimisation extérieure ne porte que sur les paramètres "shape"
$\theta = \{r_i, \log\zeta_i, \lambda_{ij}\}$ — soit ~$113 - N = 98$ paramètres.

**Loss** :
$$\mathcal{L}(\theta) = \|w_z - \Pi_\theta\, w_z\|_{H^1}^2$$
où $\Pi_\theta = \sum_i c_i^\star(\theta) \Phi_i(\theta)$ est la projection
Tikhonov-régularisée. Pas de pénalité $\theta$-niveau pour v1.

**Optimiseur** : **L-BFGS avec K restarts indépendants** (`Optim.jl`).
- Convergence par restart : `g_abstol=1e-5`, `iterations=200`, `f_reltol=1e-8`
- Restarts diversifiés par RNG seed (et optionnellement par split
  $N_\text{centered}/N_\text{pibond}$)
- Best-of-K retenu en sortie

L'argument pour les restarts plutôt qu'AdamW : le problème est lisse,
déterministe, de petite dimension. L-BFGS converge en super-linéaire ; les
restarts couvrent la non-convexité par diversification d'init.

## 4. Paramétrisation et bornes

### 4.1 Borne sur $\log\zeta$

Box-via-sigmoid en log-space :
$$\log\zeta_i = \log\zeta_\min + (\log\zeta_\max - \log\zeta_\min)\cdot\sigma(u_i),
\quad u_i \in \mathbb{R} \text{ (paramètre libre)}$$

Bornes (calées sur Ecut, justification phase A) :
- $\zeta_\min = 10^{-2}$ (gaussienne très large, pas de problème numérique)
- $\zeta_\max = \mathrm{Ecut}/4$ — typiquement ~4 pour Ecut=15, ~7.5 pour Ecut=30, ~12.5 pour Ecut=50

### 4.2 Centre $r$ (type 2)

$r \in [0.5, 5.0]$ Å, paramétrisation directe (clamping en cas de débordement
pendant l'optim). Justification : tous les centres greedy sont dans cette
plage.

### 4.3 Coefficients polynomiaux

Pas de paramétrisation contrainte. $\lambda_{ij} \in \mathbb{R}$, optim libre.

### 4.4 z-symétrie

$\alpha_z = 0$ forcé pour préserver la symétrie $A_2''$ ($n_z$ impair). Pas
de paramètre libre en $z$.

## 5. Initialisation

Pour chaque restart $k = 1, \ldots, K$, RNG seed différent. Sampling :

| Paramètre | Distribution |
|-----------|--------------|
| $u$ (latent pour $\log\zeta$) | $\mathcal{N}(0, 1.5)$ — couvre toute la box via sigmoid |
| $r$ (type 2) | uniforme dans $[0.5, 5.0]$ Å |
| $\lambda_{ij}$ (coefs des groupes polynomiaux) | valeurs cohérentes avec table 5.2 + perturbation $\mathcal{N}(0, 0.05)$ |

(Les $c_i$ sont obtenus par variable projection à la première éval, donc
pas d'init explicite.)

Init des coefs polynomiaux aux valeurs symétrie-adaptées de la table 5.2 (par
exemple $[1, 2, 1]$ pour $n_x = 4$) avec petite perturbation : donne une
orientation physiquement cohérente, la perturbation casse la symétrie
inter-restarts pour la diversité.

## 6. Hyperparamètres

| Hyperparam | Default | Notes |
|------------|---------|-------|
| $N_\text{centered}$ | 8 | varier selon expériences |
| $N_\text{pibond}$ | 7 | $N = N_\text{centered} + N_\text{pibond}$ |
| $K$ (restarts) | 3 (smoke) / 5–10 (run) | cohérent avec section 9 |
| $\varepsilon$ (Tikhonov) | $10^{-8}$ | repris de phase A |
| $\zeta_\min$ | $10^{-2}$ | |
| $\zeta_\max$ | $\mathrm{Ecut}/4$ | |
| Ecut | 15 (smoke) / 30 (local) / 50 (cluster) | |
| `max_xy_order`, `max_z_order` | 3, 3 | comme phase A |

## 7. Backend

- **Julia 1.11+**, `Optim.jl` (L-BFGS)
- **AD : Zygote** par défaut. Si trop lent (forward+backward sur FFT
  supercell), bascule vers Enzyme. Pas d'optim prématurée.
- Tikhonov inner solve : closed form `(S + εI) \ Γ`. Compatible Zygote
  (linear solve adjoint via ChainRulesCore).
- Discrétisation : conserve la base plane-wave DFTK existante
  (`PlaneWaveBasis` + `cell_to_supercell`). Le wannier est lu depuis le
  JSON, le basis_supercell est reconstruit avec mêmes params.

## 8. Outputs et persistance

Pour chaque run, un dossier `workflow/joint_outputs/run_<timestamp>/`
contient :

- `config.json` — tous les hyperparams + seed maître
- `restart_<k>.json` — pour chaque $k$ : init params, loss curve, params
  finaux, loss finale, wallclock
- `best.json` — meilleur restart
- `S_final.csv` — matrice de recouvrement de la base finale
- `summary.txt` — comparatif : best vs greedy baseline, distribution des
  $\log\zeta$, top paires recouvrantes

Diagnostics réutilisés depuis phase A :
- `cond(S)`, $\sigma_\min$, $\sigma_\max$
- Top 10 paires de recouvrement normalisé
- Histogramme des $\log\zeta$ finaux

## 9. Tests et validation

### 9.1 Smoke test (Ecut=15)

- $N=10$, $K=3$, max 50 itérations L-BFGS par restart
- Vérifie : pipeline entier tourne, AD compile, loss décroît
- Cible : loss finale **< 30%** (juste pour valider l'infra)
- Coût : <5 min sur laptop

### 9.2 Run baseline (Ecut=15)

- $N=15$, $K=5$, paramètres par défaut
- Compare directement au greedy@15 de phase A (8.6%)
- Cible : **< 8% H¹** (battre le greedy avec mêmes hyperparams discrets)
- Coût attendu : ~30 min laptop

### 9.3 Run cible (Ecut=30, local)

- $N=15$, $K=10$
- Cible : **≤ 5% H¹** (matcher Bakhta+18)
- Coût attendu : ~3h laptop

### 9.4 Run précision (Ecut=50, cluster)

- $N=15-20$, $K=10$
- Cible : descendre vers 1–2% H¹ si possible
- Hors phase B stricte mais bon à avoir pour la phase aval TB

## 10. Reproductibilité

- Seed maître en input
- Versioning git du code + hyperparams sauvegardés avec outputs
- Convention : `git rev-parse HEAD` enregistré dans `config.json`

## 11. Hors scope (v1)

- Adaptation auto de $N$ pendant l'optim (régul L1 sur $c$, ou critère
  d'arrêt sur amélioration relative)
- Stochastique / mini-batching : non pertinent (cible unique)
- AdamW ou autres optimiseurs ML (option si L-BFGS+restarts ne suffit pas)
- Variation de $N_\text{centered}/N_\text{pibond}$ comme hyperparam :
  exploration manuelle d'abord
- Pénalités $\theta$-niveau (sur $\log\zeta$, conditionnement de S, etc.)
- AD via Enzyme (option de second tour si Zygote trop lent)
- GPU : la base plane-wave est petite, pas pertinent pour le moment
- Reformulation atomique (proposition B du brainstorming) : différée

## 12. Questions ouvertes

À documenter quand on aura les premiers résultats :

- Les restarts convergent-ils tous au même optimum, ou observe-t-on des
  bassins distincts ? (réponse pratique sur la non-convexité)
- L'init aux valeurs de la table 5.2 perturbées est-elle un bon biais, ou
  un init purement aléatoire (e.g., $\lambda_{ij} \sim \mathcal{N}(0, 1)$)
  donne-t-il de meilleurs résultats ?
- Au-delà de quel $N$ le fit conjoint sature-t-il ? (analogue à la question
  greedy mais en mieux)
- L'erreur descend-elle sous le seuil de précision FFT (3% pour Ecut=15) ?

## 13. Livrables

Code :
- `src/joint_optim/` — nouveau module
  - `loss.jl` — calcul de $\mathcal{L}(\theta)$ via variable projection
  - `parametrization.jl` — encodage box-sigmoid, init, conversion
    paramètres → SAGTOs
  - `runner.jl` — boucle de restarts, sauvegarde
- `workflow/joint_compression.jl` — script d'entrée
- Tests unitaires : sigmoid encoding round-trip, gradient checking, loss
  monotone décroissante sur cas trivial

Notes :
- `notes/04-bilan-phase-B.md` — rédigée *après* les premiers runs

## 14. Références

- `notes/00-analyse-et-plan.md` — cadrage projet
- `notes/02-bilan-phase-A.md` — fondement du design phase B
- [LV1] Cancès, Dusson, Kemlin, Vidal. *On basis set optimisation in
  quantum chemistry*. ESAIM Proc. Surv. 73 (2023). → précédent fit conjoint
- [Bak+18] Bakhta et al. — méthode greedy de référence
- Golub, Pereyra. *The differentiation of pseudo-inverses and nonlinear
  least squares problems whose variables separate*. SIAM J. Numer. Anal. 10
  (1973). → fondement de variable projection
