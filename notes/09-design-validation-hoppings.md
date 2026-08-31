# Design — validation des hoppings t(R≠0)

Étape suivante de la chaîne de validation TB (notes 05 §4 et 08 §10).
Design validé au brainstorm du 29 août 2026 : **approche A puis B**.

## Objectif

Mesurer l'erreur réelle du Wannier gaussien-compressé sur les hoppings
⟨w₀|·|w_R⟩, R≠0, contre la référence plane-wave. C'est **la mesure décisive**
identifiée par la note 08 §4 :

> La superconvergence de Galerkin (erreur cinétique ≈ (erreur H¹)²) ne vaut que
> pour l'élément diagonal. Pour un hopping, les deux fonctions sont projetées
> sur des espaces différents : plus de structure variationnelle, donc erreur du
> **premier** ordre — attendue ~12 % (Ecut=50), pas ~1.5 %.

Si la prédiction se vérifie, la question §6-Q2 de la note 05 est tranchée : la
métrique H¹ est trompeusement bonne sur l'on-site et honnête sur les hoppings
→ la loss de la phase B doit être réorientée vers les hoppings (elle est
différentiable de bout en bout à travers les intégrales, donc c'est faisable).

## Contexte du brainstorm (décisions amont)

- **Livrable phare** : bandes plates ab initio à l'angle magique, TB sans
  paramètre empirique. Axe secondaire : exploiter la **différentiabilité** des
  intégrales GTO (forces de Hellmann-Feynman pour la relaxation, ∂t/∂R =
  éléments électron-phonon, sensibilités d(bande)/dθ, dérivation asymptotique
  de BM à la CGG23, loss hopping pour la phase B).
- **Potentiel en TB** : pas de SCF ; V_TBG ≈ V₁+V₂ (superposition de potentiels
  monocouche, note 08 §9), validée plus tard par un SCF grand angle.
- **Argument GTO vs Fourier** : intralayer, gain quantitatif ; **interlayer
  twisté, gain qualitatif** (les grilles G des deux couches tournées sont
  incommensurables — en Fourier il faudrait la grille moiré ~10⁸-10⁹
  coefficients ; une gaussienne se tourne analytiquement, coût indépendant de
  l'angle). Plus la différentiabilité, impossible en Slater-Koster/Fourier pur.

## Approche A — overlap + cinétique, translation rigide

Aucune brique nouvelle : étendre `workflow/validate_tb_onsite.jl` en
`workflow/validate_tb_hoppings.jl`.

### Les deux chemins

| Côté | Translation | Intégrale |
|---|---|---|
| Référence (plane-wave) | phase e^{−iG·R} sur les coefficients Fourier | somme sur grille (machinerie `_check_L²_dot_precision` / `compare_onsite_kinetic`) |
| Gaussien (compressé) | `translate(Wc, R)` (`src/CompressedWannier.jl:33`) | `integral(Wc₁, Wc₂; type=:overlap/:kinetic)` (`src/integrals/integrals.jl:63`), natif GaIn-free |

### Échantillonnage des R (validé)

- **Intralayer** (3 hoppings) : a₁, a₁+a₂, 2a₁ — premiers, seconds, troisièmes
  voisins de cellule.
- **Pseudo-interlayer** (3-4 hoppings) : translation rigide d = 3.35 Å en z,
  avec décalage latéral δ ∈ {0 (AA-like), (a₁+a₂)/3 (AB-like), a₁/2
  (intermédiaire)}. Pas de rotation à ce stade — la translation z suffit pour
  sonder la queue p_z hors-plan (drapeau note 05 §5).

### Quantités mesurées

Pour chaque R et chaque compressé (`data/compressed_wannier_55_H1.json`
Ecut=15, `data/compressed_wannier_55_H1_Ecut50.json` Ecut=50) :

- S(R) = ⟨w₀|w_R⟩ et T(R) = ½⟨∇w₀,∇w_R⟩, référence vs gaussien ;
- erreur relative |ΔX|/|X_ref| par hopping ; attention aux petits
  dénominateurs (hoppings lointains) → rapporter aussi l'erreur absolue
  normalisée par T(0) ;
- tableau récapitulatif erreur vs distance, comparé à l'on-site (0.41 % /
  1.89 %) et à la prédiction premier ordre (≈ erreur H¹ : 5.59 % / 12.38 %).

### Critères de succès

L'expérience est concluante si elle **discrimine** les deux régimes :
- erreur hopping ≫ erreur on-site et ~erreur H¹ → prédiction confirmée,
  réorientation de la loss phase B actée ;
- erreur hopping ~ erreur on-site → superconvergence plus robuste que prévu,
  le fit H¹ suffit, la chaîne TB continue telle quelle.

Dans les deux cas on gagne ; le danger serait un résultat intermédiaire brouillé
par les petits dénominateurs — d'où la double métrique relative/absolue.

## Approche B — hopping physique complet (après A)

Ajouter le terme potentiel pour t(R) = ⟨w₀|−½∇²+V_KS|w_R⟩ en eV, comparable à
la littérature (t₁ ≈ −2.7 eV premier voisin du graphène).

- Câbler un `scfres` monocouche cohérent cellule↔supercellule (L-2b, note 08
  §10) ; V_KS via `potential_scalar_prod` (déjà GaIn-free, Fourier).
- Mêmes R que A ; référence = même intégrale sur le Wannier plane-wave.
- Ancrage absolu : |t₁^{gto} − t₁^{litt}| et le ratio t₂/t₁ (sensible à la
  portée de la queue).

B est le prérequis de T-band (bandes monocouche) et de toute la suite ; son
design fin (obtention du scfres, conventions de cellule) sera précisé dans le
plan une fois A conclu.

## Hors scope

- Rotation des gaussiennes (interlayer twisté réel) — étape TBG, après T-band.
- Nouvelle compression (greedy ζ_max, phase B Tasks 10-12) — thread parallèle ;
  on valide d'abord avec les compressés de février.
- Terme deux corps, relaxation, électron-phonon — perspectives, notées §Contexte.
