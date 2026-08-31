# Design — loss hopping-ciblée (réorientation phase B)

Suite directe du bilan [11](11-bilan-validation-hoppings.md) et de la
discussion critique du 31 août 2026. Décision amont validée : la compensation
cinétique/potentiel éventuelle dans t(R) complet ne serait pas transportable au
TBG (V₁+V₂ n'est le potentiel propre d'aucune des deux Wanniers), donc **T et S
doivent être bons en soi** avant tout Hamiltonien avec potentiel. L'approche B
(t(R) complet, scfres L-2b) est repoussée *après* ce travail, en validation
finale ancrée littérature.

## 1. Diagnostic et mécanisme

Le greedy H¹ (phase A) place son erreur là où la norme H¹ paie : le cœur de la
Wannier. Or S(R) et T(R) ne sondent que la région de recouvrement entre
voisins (2–3.6 Å), quasi invisible dans H¹. Résultat (note 11) : erreur
absolue ~η·T(0) sur des hoppings minuscules → ordre 1 en relatif, signes faux,
orthonormalité MLWF violée (S_GTO(a₁) ≈ 7e-3 vs 0 exact).

Identifiabilité : {S(R), T(R)} même pour *tous* les R ne déterminent que
l'autocorrélation de w (module de ŵ, pas sa phase). Les termes hoppings ne
peuvent donc pas *remplacer* l'approximation de la fonction — ils
**redistribuent spatialement** l'erreur d'une approximation qui reste
l'ossature de l'objectif. La notion de base minimale (courbe erreur vs N)
demeure.

## 2. Loss

$$\mathcal{L}(p) = \underbrace{\|w - w_{\text{GTO}}(p)\|_{H^1}^2}_{\text{forme / identifiabilité}}
+ \mu \sum_{R \in \mathcal{T}} \omega_R \left[ \Delta T(R)^2 + \Delta S(R)^2 \right]
+ \nu \sum_{R_{\text{intra}} \in \mathcal{T}} S_{\text{GTO}}(R)^2$$

avec ΔX(R) = X_GTO(R) − X_ref(R), valeurs norm-corrected des deux côtés
(convention `compare_hopping`).

- **ω_R = 1/T(0)²** : adimensionne, aligne les termes sur la métrique
  err/T(0) du bilan 11.
- **μ** : scan log sur 3-4 décades, smoke Ecut 15, sélection par les critères
  §5 sur le jeu de *validation*.
- **ν** : initialisé à μ, augmenté si |S_GTO(a₁)| plafonne au-dessus du
  critère (3).
- **v2 documentées, non implémentées** : (i) lagrangien augmenté si la
  pénalité ν plafonne ; (ii) formulation riemannienne de la
  "shift-orthonormalité" — la contrainte S_GTO(R_intra; p) = 0 n'est pas une
  contrainte de Stiefel sur p (égalités non linéaires, colonnes = translatées
  d'une même fonction à paramètres partagés, la structure produit ne factorise
  pas), mais une formulation propre est peut-être possible : **décision
  déléguée à Laurent** (partie 1 de sa thèse : optim sur Stiefel/Flag,
  bibliothèques R-LBFGS compatibles ForwardDiff connues de lui). Morceau
  méthodologique potentiel en soi.

## 3. Véhicule : option (c), hybride

Init greedy H¹ (phase A, inchangé) → raffinement joint-optim (phase B,
L-BFGS + ForwardDiff) avec la loss augmentée. Ajouter les termes hoppings à
`joint_loss` est architecturalement trivial : les intégrales gaussiennes
S_GTO(R), T_GTO(R) sont analytiques et différentiables (machinerie
`compare_hopping` / `integral` + `translate`).

**Prérequis absorbés dans le périmètre** : Tasks 10 (K-restarts +
persistance) et 11 (smoke workflow) de la phase B, en pause depuis mai
(cf. notes 03/04 et ancien RESUME-HERE dans l'historique git). En cas d'échec
de la réorientation, cette machinerie sert telle quelle à l'approche B.

**Coût** : les termes ajoutés sont marginaux devant le terme H¹ sur grille
(goulot ForwardDiff existant : ~7.5 s/iter à 20 params, ~42 s/iter extrapolé à
113, cf. benchmark Task 9). Budgets de mai inchangés en ordre de grandeur.

## 4. Jeux de R et cibles

- **Entraînement 𝒯** : jeu note 09 moins 2a₁ (contaminé par les images
  périodiques de la supercellule 5×5) : a₁, a₁+a₂, inter AA (z=3.35 Å),
  inter AB-like (z + (a₁+a₂)/3), inter mid (z + a₁/2). Terme ν sur a₁ et
  a₁+a₂.
- **Validation** (jamais dans la loss) : cœur déterministe
  {z + a₁/4, z'=3.0 Å (δ=0), z''=3.7 Å (δ=0), a₂} + ~5 R aléatoires dans la
  couronne 2–4 Å, `Random.seed!(SEED)` avec `SEED` constante nommée du script
  (changer de seed explicitement = test de robustesse du verdict).
- **Cibles précalculées** : script dédié calculant S_ref(R), T_ref(R) (côté
  plane-wave, phases e^{−iG·R}) une fois pour toutes → JSON dans `data/`
  (par Ecut). Zéro plane-wave dans la boucle d'optim pour les termes hoppings.
- Ecut 15 = smoke (mise en place, hard fails) ; Ecut 50 = runs de sélection et
  verdict. Rappel note 11 : la référence T_ref n'est pas convergée en Ecut
  (facteur 3 sur T_ref(a₁) entre 15 et 50) — les cibles sont *par Ecut*, on ne
  compare jamais des runs à Ecut différents.

## 5. Critères de succès

Mesurés sur le jeu de **validation**, Ecut 50, chacun **loggé indépendamment**
(valeur + PASS/FAIL par critère — voir lequel mord en premier ; si (4) casse
seul, le goulot est le budget N, pas la loss) :

1. relerr T(a₁) ≤ **5 %** (actuellement 69 %) — retour au premier ordre ~η ;
2. **signes corrects** sur tous les T de validation avec |T_ref| > 5e-4 Ha
   (en dessous : bruit de la référence) ;
3. |S_GTO(a₁)| ≤ **1e-3** (actuellement 7.2e-3 à Ecut 15 / 1.2e-3 à Ecut 50) ;
4. anti-régression : erreur H¹ ≤ **1.5×** celle du run jumeau μ=0.

**Comparaison contrôlée obligatoire** : chaque run va par paire (μ=0 vs μ>0),
même init greedy, même N, même budget d'itérations. C'est la mesure qui valide
ou invalide la réorientation — pas les valeurs absolues seules.

Échec des critères (1-3) à μ optimal et ν poussé → le budget N=10-15 fonctions
est insuffisant pour la queue → rouvrir le design compression (plus de
fonctions, ζ_min plus petits) *avec* cette loss comme objectif — pas avant.

## 6. Positionnement bibliographique

Recherche du 31 août 2026 (agent web). Verdict : **combinaison nouvelle
d'ingrédients standards** ; le gap précis — orbitales localisées *analytiques*
fittées pour que leurs S(R), T(R) induits collent à une référence ab initio —
semble inoccupé (le downfolding TBG fitte des tables/formes fonctionnelles de
t(R), pas des orbitales). Deux confirmations indépendantes du diagnostic
"norme ≠ observables" : ABACUS (spillage L² → ajout d'un terme gradient) et la
lignée Dunning/bases property-targeted.

Références prioritaires :
1. Bakhta, Cancès, Cazeaux, Fang, Kaxiras 2018 (arXiv:1712.02996) — baseline
   Wannier→GTO, loss norme Sobolev pure (notre point de départ) ;
2. Zhou, Gull, Zgid 2021 (arXiv:2105.00589) — analogue le plus proche : GTO
   optimisés contre référence plane-wave, loss mixte énergie + bandes +
   conditionnement de S ;
3. Fang & Kaxiras 2016 (PRB 93, 235153) — formes fonctionnelles t(r,θ) pour
   bilayer, la référence de ce que nous *ne* faisons *pas* (tables) ;
4. Koshino et al. 2018 (PRX 8, 031087) — downfolding Wannier standard du TBG ;
5. Lin, Ren, He 2021 (ABACUS, generalized spillage) — l'analogue NAO le plus
   proche côté loss.

(Laurent relit ces articles ; le positionnement pourra être affiné dans la
note de bilan.)

## 7. Hors scope

- Rotation interlayer réelle (étape TBG) ;
- Approche B — t(R) complet avec V_KS, scfres L-2b (piste :
  `TwistedBilayerGraphene.jl`, code de thèse de Laurent calculant le potentiel
  monolayer) — **après** ce travail, comme validation finale ;
- Re-design du greedy phase A ;
- Task 12 de l'ancien plan phase B (comparaison baseline H¹) — remplacée par
  la comparaison contrôlée μ=0/μ>0 du §5.
