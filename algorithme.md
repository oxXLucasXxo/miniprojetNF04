# Algorithme d’intégration par la méthode de Monte Carlo

## 1. Principe

L’algorithme estime l’intégrale d’une fonction polynomiale sur un intervalle `[a, b]` en utilisant un tirage aléatoire de points. Il détermine d’abord un rectangle qui englobe la courbe, puis tire un grand nombre `N` de points dans ce rectangle pour évaluer l’aire sous la courbe.

## 2. Fonctions auxiliaires

### FONCTION EvaluerPolynome(coefficients, x)
*   **Rôle** : Calcule la valeur du polynôme `P(x)` pour une valeur `x` donnée.
*   **Entrées** :
    *   `coefficients` : Liste des coefficients du polynôme, de `a_0` à `a_n`.
    *   `x` : Point d’évaluation (réel).
*   **Sortie** : La valeur `P(x)` (réel).
*   **Pseudo-code** :
    ```
    resultat ← 0
    POUR i DE 0 A degre(coefficients) FAIRE
        resultat ← resultat + coefficients[i] * (x PUISSANCE i)
    FIN POUR
    RETOURNER resultat
    ```

### FONCTION CalculerDerivee(coefficients)
*   **Rôle** : Calcule les coefficients du polynôme dérivé `P'(x)`.
*   **Entrée** :
    *   `coefficients` : Liste des coefficients de `P(x)`.
*   **Sortie** : Liste des coefficients de `P'(x)`.
*   **Pseudo-code** :
    ```
    SI degre(coefficients) < 1 ALORS
        RETOURNER [0]
    FIN SI

    coeffs_derivee ← NOUVELLE LISTE
    POUR i DE 1 A degre(coefficients) FAIRE
        AJOUTER (coefficients[i] * i) à coeffs_derivee
    FIN POUR
    RETOURNER coeffs_derivee
    ```

### FONCTION TrouverRacines(coefficients)
*   **Rôle** : Trouve les racines réelles d'un polynôme. La résolution analytique étant complexe pour un degré élevé, on suppose l'existence d'une fonction numérique qui effectue cette tâche.
*   **Entrée** :
    *   `coefficients` : Liste des coefficients du polynôme.
*   **Sortie** : Liste des racines réelles.

## 3. Algorithme principal

### FONCTION EstimerIntegraleMonteCarlo(a, b, coeffs_polynome, N)
*   **Rôle** : Estime l'intégrale de `P(x)` sur `[a,b]` avec `N` points.
*   **Entrées** :
    *   `a`, `b` : Bornes de l'intervalle (réels).
    *   `coeffs_polynome` : Coefficients du polynôme.
    *   `N` : Nombre de points pour le tirage (entier > 0).
*   **Sortie** : Valeur estimée de l'intégrale (réel).

*   **Pseudo-code** :

    ```
    // --- 1. Déterminer le rectangle englobant ---

    // a. Trouver les points où la dérivée s'annule
    coeffs_derivee ← CalculerDerivee(coeffs_polynome)
    racines_derivee ← TrouverRacines(coeffs_derivee)

    // b. Identifier les points critiques pour les extrema sur [a,b]
    points_critiques ← NOUVELLE LISTE avec {a, b}
    POUR chaque racine DANS racines_derivee FAIRE
        SI a ≤ racine ≤ b ALORS
            AJOUTER racine à points_critiques
        FIN SI
    FIN POUR

    // c. Évaluer le polynôme aux points critiques pour trouver les extrema locaux
    valeurs_extremes ← NOUVELLE LISTE
    POUR chaque point DANS points_critiques FAIRE
        AJOUTER EvaluerPolynome(coeffs_polynome, point) à valeurs_extremes
    FIN POUR

    // d. Déterminer les bornes verticales du rectangle
    y_min_f ← MINIMUM(valeurs_extremes)
    y_max_f ← MAXIMUM(valeurs_extremes)
    y_min_rect ← MINIMUM(0, y_min_f)
    y_max_rect ← MAXIMUM(0, y_max_f)

    aire_rectangle ← (b - a) * (y_max_rect - y_min_rect)
    SI aire_rectangle = 0 ALORS
        RETOURNER 0 // Cas où f(x)=0 sur [a,b]
    FIN SI

    // --- 2. Tirage de Monte Carlo ---

    somme_contributions ← 0
    POUR i DE 1 A N FAIRE
        // Tirage aléatoire uniforme dans le rectangle
        x_aleatoire ← ALEATOIRE_REEL(a, b)
        y_aleatoire ← ALEATOIRE_REEL(y_min_rect, y_max_rect)

        valeur_f ← EvaluerPolynome(coeffs_polynome, x_aleatoire)

        // Calcul de la contribution du point
        SI y_aleatoire ≥ 0 ET y_aleatoire ≤ valeur_f ALORS
            somme_contributions ← somme_contributions + 1
        SINON SI y_aleatoire < 0 ET y_aleatoire ≥ valeur_f ALORS
            somme_contributions ← somme_contributions - 1
        FIN SI
    FIN POUR

    // --- 3. Calcul final ---

    integrale_estimee ← (somme_contributions / N) * aire_rectangle
    RETOURNER integrale_estimee
    ```

## 4. Analyse de la complexité

Soit `n` le degré du polynôme et `N` le nombre de points tirés.

*   **`EvaluerPolynome`** : Complexité en **O(n)**. Un calcul naïf de puissance peut être plus lent, mais une méthode comme celle de Horner l'optimise.
*   **`CalculerDerivee`** : Complexité en **O(n)**.
*   **`TrouverRacines`** : La complexité dépend de la méthode numérique utilisée. C'est une étape coûteuse mais exécutée une seule fois. Notons sa complexité `C_racines(n-1)`.
*   **Détermination du rectangle** :
    *   Calcul des points critiques : `O(n)` (dérivée) + `C_racines(n-1)`.
    *   Évaluation aux points critiques : Il y a au plus `n+1` points critiques (`a`, `b` et `n-1` racines). L'évaluation de chacun prend `O(n)`. Total : `(n+1) * O(n) = O(n^2)`.
    *   Complexité de cette étape : **O(n^2 + C_racines(n-1))**.
*   **Tirage de Monte Carlo** :
    *   La boucle principale s'exécute `N` fois.
    *   À chaque itération, l'appel à `EvaluerPolynome` prend `O(n)`.
    *   Complexité de cette étape : **O(N * n)**.

**Complexité totale de l'algorithme** : `O(n^2 + C_racines(n-1) + N * n)`.

Pour un `N` très grand (ce qui est le cas d'usage de la méthode), le terme `O(N * n)` devient dominant. La complexité est donc **pratiquement linéaire par rapport au nombre de points `N` et au degré du polynôme `n`**.
