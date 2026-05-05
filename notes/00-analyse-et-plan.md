# Wannier2GTO — analyse et plan de relance

Note de cadrage rédigée à la reprise du projet (mai 2026, après pause depuis 2024).

## 1. Contexte rappelé

Pour pouvoir faire du tight-binding sur graphène bicouche twisté (TBG) en
contournant le coût des intégrales bi-électroniques sur les fonctions de Wannier
brutes, on projette la Wannier $p_z$ du graphène monocouche sur une base finie
$\Phi^{\mathrm{SA}} = \{\phi_i^{\mathrm{SA}}\}_{i \le n}$ de SAGTOs (Symmetry-Adapted
Gaussian-Type Orbitals), $A_2''$-symétriques sous $D_{3h}$. Les ERI peuvent alors
être calculées analytiquement via GaIn (lib C++ d'I. Duchemin).

Méthode : compression greedy de Bakhta–Cancès–Cazeaux–Fang–Kaxiras 2018, adaptée
à $D_3$ pour la $p_z$ du graphène. À l'itération $n$ :

$$\phi_n^{\mathrm{SA}} \in \arg\min_{\phi^{\mathrm{SA}} \in \mathcal{B}} \|r_n - \phi^{\mathrm{SA}}\|_{H^s}^2,
\qquad C_n = S_n^{-1} X_n
\quad\text{(coefficients de projection } H^s\text{-orthogonale)}.$$

Détails dans le manuscrit, section 5.4.

## 2. Diagnostic du code actuel

Architecture (`src/Wannier2GTO.jl` + sous-modules) :

- `compress_graphene_pz_wannier` : boucle greedy, `Optim.jl` + `ConjugateGradient` + `ForwardDiff`
- Variables d'optim **extérieure** : `log ζ` (1 paramètre) si centré, ou `[center, log ζ]` (4 paramètres) sinon
- Niveau **intérieur** : combinaison linéaire optimale des SAGTOs pour ordres polynomiaux fixés, calculée en forme fermée par `S \ Γ` (overlap matrix \ projection sur le résidu)
- Itérations : ajout d'une SAGTO par tour, alternance `D3_sym` / non-`D3_sym`

Forces :
- Symétries explicitement encodées (table 5.2 : polynômes $D_3$-invariants)
- Closed-form pour les coefficients à chaque itération (chaque sous-problème reste petit)

Limitations identifiées (déjà notées en thèse 5.4.2.3) :
- **Bloquant : conditionnement de $S$** — empêche d'augmenter la base au-delà
  de ~10 fonctions (155 paramètres totaux), résiduel $H^1$ plateau à ≈12%
- AD en mode forward (lent à plus haute dim)
- Pas de procédure de cure du conditionnement (TODO ligne 119 de `compression.jl`)

## 3. Réponses aux questions de cadrage

### 3.1. Le code s'inscrit-il dans un cadre d'optim ML ?

**Partiellement, dans la formulation actuelle** :

- ✅ AD utilisée, problème non-convexe
- ❌ **1–4 paramètres** par itération extérieure — très loin du régime ML
- ❌ Pas de stochasticité naturelle ; sous-problème intérieur quadratique fermé
- ❌ Le bloquant *réel* est numérique (cond$(S)$), pas algorithmique

→ Migrer vers Adam/AdamW dans cette formulation **ne résoudra pas** le problème principal.

**Reformulation B (alternative)** : abandonner le greedy, fitter conjointement
tous les paramètres (centres, spreads, coefficients) sur
$\mathcal{L}(\theta) = \|w_z - \sum_i c_i \phi_i^{\mathrm{SA}}(\theta)\|_{H^s}^2 + \mathcal{R}(\theta)$.
Là on a ~150–500 paramètres → ML pertinent (Flux + AdamW + scheduler).
C'est le paradigme « trainable basis set » (cf. [LV1] Cancès–Dusson–Kemlin–Vidal 2023).

### 3.2. PINS / Herglotz-NET donne-t-il un meilleur cadre ?

**Non, pas directement.** Différences fondamentales :

| | Herglotz-NET | Wannier2GTO |
|---|---|---|
| Domaine | $\mathbb{S}^2$ | $\mathbb{R}^3$ avec périodicité 2D |
| Base | harmoniques sphériques (via PE Herglotz) | SAGTOs $A_2''$ centrées |
| Représentation | INR (MLP continu) | combinaison de gaussiennes-polynômes |
| Symétrie encodée | $\mathrm{SO}(3)$ via PE | $D_{3h}$ via construction |

L'analogie « PE adapté à la symétrie » existe mais reste superficielle. Surtout,
remplacer la base par un MLP **casse le point fort du projet** :
l'analyticité des intégrales ERI gaussiennes, qui est la raison d'être de la
compression. Donc on n'adapte pas Herglotz-NET ici.

### 3.3. Refactor avec Flux.jl ?

Pertinent **uniquement dans la voie B** (fit conjoint). Pour la voie A
(diagnostic conditionnement), `Optim.jl` + AD reverse suffit largement.

## 4. Plan retenu : A puis B

### Phase A — Réparer le conditionnement (priorité immédiate)

Voir `01-plan-A-conditionnement.md`. Esprit : 10 lignes de diagnostic,
puis tester 2–3 remédiations standard (Löwdin symétrique, Cholesky pivotée,
Tikhonov). Pas de Flux ici. Critère de succès : pouvoir grandir la base
au-delà de 15–20 fonctions sans explosion de cond$(S)$.

### Phase B — Reformulation conjointe ML (à brainstormer après A)

Sera cadrée *après* A, en fonction de ce qu'on apprend (notamment :
quels termes tirent vers le mauvais conditionnement, sont-ils
indépendamment éliminables ou inhérents ?). Brainstorming dédié
en temps voulu. Pistes pré-identifiées (à challenger) :

- Loss conjoint $\|w_z - \sum_i c_i \phi_i\|_{H^1}^2$ + régularisation
  contre la dégénérescence (pénaliser cond$(S)$ ou ajouter Tikhonov)
- Initialisation depuis un greedy court (warm-start pour échapper aux mauvais minima)
- Reverse-mode AD (Zygote ou Enzyme), AdamW + scheduler
- Garder l'invariance $D_3$ par construction des SAGTOs (pas via régul)

## 5. Références clés

- **[Bak+18]** Bakhta, Cancès, Cazeaux, Fang, Kaxiras. *Compression of Wannier
  functions into Gaussian-type orbitals*. Comp. Phys. Comm. 230 (2018), 27–37.
  → méthode greedy originelle.
- **Manuscrit Vidal 2024**, section 5.4 (p. 158–164) : adaptation $D_3$ pour
  graphène, résultats préliminaires, perspectives non finalisées.
- **[LV1]** Cancès, Dusson, Kemlin, Vidal. *On basis set optimisation in quantum
  chemistry*. ESAIM Proc. Surv. 73 (2023), 107–129. → précédent de fit conjoint.
- **[Leh19b]** Lehtola. *Curing basis set overcompleteness with pivoted Cholesky
  decompositions*. arXiv:1911.10372. → remède standard pour A.
- **[Duc]** Duchemin. *GaIn — Gaussian Integral library*. → cible aval.
- **Herglotz-NET** (Hanon et al. 2026, arXiv:2502.13777) : non retenu pour ce
  projet (cadre incompatible).
