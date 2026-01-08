# Intégration numérique d’une fonction polynomiale (méthode de Monte Carlo)

## Objectif

Implémenter en Python une méthode **approximative** de calcul de l’intégrale d’une fonction réelle, en se limitant au cas des **fonctions polynomiales de degré inférieur ou égal à 20**, définies sur un intervalle réel ([a,b]).

La méthode utilisée repose sur un **tirage aléatoire uniforme de points** dans un rectangle englobant la courbe de la fonction, afin d’estimer l’aire algébrique sous la courbe.

---

## Fonction étudiée

La fonction considérée est un polynôme réel :

[
f(x) = a_0 + a_1 x + a_2 x^2 + \dots + a_n x^n
]

avec :

* (n \leq 20)
* coefficients réels
* intervalle d’intégration ([a,b]), avec (a < b)

---

## Rectangle d’intégration

Le rectangle utilisé pour le tirage aléatoire doit **englober entièrement la courbe** de (f) sur ([a,b]).

Il est défini par :

* l’intervalle horizontal : ([a,b])
* l’intervalle vertical :
  [
  \left[\min(0, \min_{x \in [a,b]} f(x)),; \max(0, \max_{x \in [a,b]} f(x))\right]
  ]

La détermination des **extrema du polynôme sur ([a,b])** est donc une étape nécessaire.

---

## Méthode de Monte Carlo

1. On tire aléatoirement des points ((x,y)) selon une **loi uniforme** dans le rectangle.
2. Chaque point contribue :

   * **+1** si (0 \le y \le f(x))
   * **−1** si (f(x) \le y \le 0)
   * **0** sinon
3. L’intégrale est estimée par :

[
\int_a^b f(x),dx \approx
\left(\frac{\sum \text{contributions}}{N}\right)
\times \text{aire du rectangle}
]

où (N) est le nombre de points tirés.

Le nombre de tirages est **choisi par l’utilisateur** et peut être très grand.

---

## Données d’entrée

L’utilisateur fournit :

1. Les bornes de l’intervalle ([a,b])
2. Les coefficients du polynôme
3. Le nombre de points aléatoires (N) pour l’estimation

---

## Résultats attendus

Le programme doit :

* calculer les extrema du polynôme sur ([a,b])
* afficher ou tracer :

  * la courbe du polynôme
  * la surface correspondant à l’intégration
* estimer numériquement l’intégrale par la méthode de Monte Carlo
* afficher la valeur approchée de l’intégrale

---

## Contraintes et remarques

* Le langage utilisé est **Python**
* L’algorithme doit être générique (tout polynôme de degré ≤ 20)
* La précision dépend du nombre de tirages
* Les cas où la fonction change de signe doivent être correctement pris en compte
