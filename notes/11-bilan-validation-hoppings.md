# Bilan — validation hoppings, approche A

Exécution du plan [10](10-plan-validation-hoppings-A.md), design
[09](09-design-validation-hoppings.md). Logs :
`workflow/diag_outputs/hoppings_Ecut{15,50}.log`. Valeurs norm-corrected
(⟨·⟩/‖w‖² de part et d'autre).

## 1. Résultats

### Ecut = 15 (H¹ = 5.59 %)

| R | \|R\| (Å) | S_ref | S_gto | T_ref | T_gto | relerr_T | err/T(0) |
|---|---|---|---|---|---|---|---|
| 0 (check) | 0 | 1 | 1 | 1.4721 | 1.4814 | **0.63 %** | 0.63 % |
| intra a₁, a₁+a₂ | 2.64 | ~0 (2e-16) | 7.2e-3 | 1.301e-2 | 2.047e-2 | **57 %** | 5.1e-3 |
| intra 2a₁ | 5.28 | ~0 | −1.7e-3 | −7.8e-5 | 2.9e-4 | 479 % | 2.5e-4 |
| inter AA | 3.35 | −4.95e-2 | −4.82e-2 (2.5 %) | −8.9e-4 | **+1.3e-3** | 244 % | 1.5e-3 |
| inter AB-like | 3.46 | −3.76e-2 | −3.92e-2 (4.4 %) | −2.1e-5 | −1.2e-3 | 5930 % | 8.3e-4 |
| inter mid | 3.60 | −2.59e-2 | −2.95e-2 (14 %) | 1.0e-3 | **−2.1e-3** | 306 % | 2.2e-3 |

### Ecut = 50 (H¹ = 12.38 %)

| R | \|R\| (Å) | S_ref | S_gto | T_ref | T_gto | relerr_T | err/T(0) |
|---|---|---|---|---|---|---|---|
| 0 (check) | 0 | 1 | 1 | 1.5051 | 1.4833 | **1.45 %** | 1.45 % |
| intra a₁, a₁+a₂ | 2.64 | ~0 | 1.2e-3 | 3.886e-2 | 1.210e-2 | **69 %** | 1.8e-2 |
| intra 2a₁ | 5.28 | ~0 | 2.5e-3 | 1.4e-4 | 6.6e-4 | 388 % | 3.5e-4 |
| inter AA | 3.35 | −5.52e-2 | −5.32e-2 (3.7 %) | −2.4e-3 | **+2.3e-4** | 109 % | 1.8e-3 |
| inter AB-like | 3.46 | −3.99e-2 | −3.85e-2 (3.5 %) | −7.4e-4 | **+1.0e-3** | 242 % | 1.2e-3 |
| inter mid | 3.60 | −2.58e-2 | −2.60e-2 (0.9 %) | 9.3e-4 | 1.0e-3 | 11 % | 6.7e-5 |

Contrôles : R=0 reproduit la note 08 (valeurs normalisées 0.63 % / 1.45 % —
la note 08 §T(0) citait 0.41 % / 1.89 % en brut ; T(0)_ref = 1.4721265 et
1.5051430 identiques à tous les chiffres) ; a₁ et a₁+a₂ strictement dégénérés
(symétrie C₃ respectée par le compressé D3-symétrisé) ; S_ref intra ~1e-16 =
orthonormalité exacte du set de Wannier.

## 2. Verdict sur la prédiction (note 08 §4)

**Branche 1 du design confirmée, et au-delà.** L'erreur hopping n'est pas
seulement du premier ordre (~H¹) : sur T(R) elle est **d'ordre 1 en relatif**
(57–69 % intra, signes faux en interlayer). La raison est un effet de
compensation : T(R≠0) est une petite différence de grandes quantités
(T(R)/T(0) ~ 10⁻²), donc une erreur absolue "premier ordre" ~ η·T(0) suffit à
tout emporter. C'est exactement ce que montre la colonne err/T(0) : 0.5–1.8 %,
soit l'ordre η·(échelle), cohérent avec la perte de superconvergence — mais
rapportée à un dénominateur minuscule.

Trois faits saillants au-delà de la prédiction :

1. **Signes faux** : en interlayer, T_gto a le mauvais signe dans 3 cas sur 6.
   Un TB construit là-dessus aurait des hoppings cinétiques qualitativement
   faux — pas juste imprécis.
2. **Violation d'orthonormalité** : S_ref(a₁) = 0 exactement (orthonormalité
   MLWF), S_gto(a₁) = 7.2e-3 (Ecut 15) / 1.2e-3 (Ecut 50). Dans un TB
   non-orthogonal, la matrice S du compressé est fausse là où la vraie est
   l'identité (entre images de réseau d'une même couche).
3. **La référence elle-même n'est pas convergée** : T_ref(a₁) = 1.30e-2
   (Ecut 15) vs 3.89e-2 (Ecut 50) — un facteur 3. Les deux wannierizations
   sont indépendantes (gauges potentiellement différentes) et le hopping
   cinétique nu est une quantité fine. Toute comparaison à la littérature
   devra se faire sur le hopping *complet* t(R) (approche B), pas sur T seul.

## 3. Mise à l'échelle physique

t₁ (premier voisin, graphène) ≈ −2.7 eV ≈ −0.1 Ha. L'erreur absolue sur le
terme cinétique intra à Ecut 50 est |ΔT(a₁)| = 2.7e-2 Ha ≈ **0.73 eV, soit
~27 % de t₁** si rien ne la compense dans le terme potentiel. La chaîne TB en
l'état ne produit pas des hoppings fiables. Réserve : H·w a une structure
(w quasi-combinaison d'états propres), une compensation partielle
cinétique/potentiel dans t(R) complet n'est pas exclue — c'est précisément ce
que l'approche B mesurera avant toute conclusion définitive.

## 4. Réserves méthodologiques

- **Images périodiques** : référence 5×5 périodique vs GTO espace libre. Pour
  2a₁ (5.28 Å vs 13.2 Å de supercellule) le plancher d'images n'est pas
  négligeable — et S_gto(2a₁) > S_gto(a₁) à Ecut 50 (2.5e-3 vs 1.2e-3) sent
  la queue de gaussienne large, à recouper. Les lignes 2a₁ sont indicatives.
- **relerr_S intra** (~1e13) : dénominateur exactement nul par orthonormalité,
  colonne sans signification là — err/T(0) et la valeur brute S_gto font foi
  (anticipé au design).
- Les hoppings cinétiques nus interlayer sont ~1e-3 Ha, à la limite du bruit
  de compression : les relerr_T interlayer (11 %–5930 %) sont peu robustes
  individuellement ; c'est le motif d'ensemble (err abs ~η·T(0), signes
  instables) qui est significatif.

## 5. Décisions

1. **Réorientation de la loss phase B : actée** (question §6-Q2 de la note 05
   tranchée). La métrique H¹ est aveugle à ce que le TB consomme. La loss doit
   cibler les hoppings — directement les quantités ⟨w₀|·|w_R⟩ sur un jeu de R
   physiques (elles sont différentiables à travers les intégrales analytiques),
   avec l'orthonormalité S(R_intra) = 0 en contrainte ou pénalité.
2. **Approche B (t(R) complet avec V_KS) reste le juge de paix** avant de
   toucher à la compression : elle dira si la compensation cinétique/potentiel
   sauve les hoppings complets, et donnera l'ancrage littérature (t₁, t₂/t₁).
   Prérequis inchangé : scfres cohérent cellule↔supercellule (L-2b, note 08
   §10).
3. La comparaison littérature sur T seul est abandonnée (référence non
   convergée en Ecut, §2.3).
