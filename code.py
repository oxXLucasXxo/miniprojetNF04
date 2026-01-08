# -*- coding: utf-8 -*-
"""
Implmentation de l'intgration numrique par la mthode de Monte Carlo
pour des fonctions polynomiales.
"""

import random
import numpy as np
from typing import List

# --- Fonctions auxiliaires ---

def evaluate_polynomial(coefficients: List[float], x: float) -> float:
    """
    Calcule la valeur du polynme P(x) pour une valeur x donne.

    Args:
        coefficients: Liste des coefficients du polynme, de a_0   a_n.
        x: Point dvaluation.

    Returns:
        La valeur P(x).
    """
    # np.polyval attend les coefficients du plus haut degr au plus bas,
    # il faut donc inverser la liste.
    return np.polyval(coefficients[::-1], x)

def calculate_derivative_coeffs(coefficients: List[float]) -> List[float]:
    """
    Calcule les coefficients du polynme driv P'(x).

    Args:
        coefficients: Liste des coefficients de P(x), de a_0   a_n.

    Returns:
        Liste des coefficients de P'(x), de a'_0   a'_{n-1}.
    """
    if len(coefficients) < 2:
        return [0.0]

    # np.polyder attend les coefficients du plus haut degr au plus bas.
    poly_obj = np.poly1d(coefficients[::-1])
    derivative_poly_obj = np.polyder(poly_obj)
    
    # On retourne les coefficients dans l'ordre a'_0, a'_1, ...
    return list(derivative_poly_obj.coeffs)[::-1]

def find_real_roots(coefficients: List[float]) -> List[float]:
    """
    Trouve les racines relles d'un polynme.

    Args:
        coefficients: Liste des coefficients du polynme, de a_0   a_n.

    Returns:
        Liste des racines relles.
    """
    if not coefficients or all(c == 0 for c in coefficients):
        return []
        
    # np.roots attend les coefficients du plus haut degr au plus bas.
    roots = np.roots(coefficients[::-1])
    
    # On ne garde que les racines relles.
    return [r.real for r in roots if np.isreal(r)]

# --- Algorithme principal ---

def estimate_integral_monte_carlo(
    a: float, 
    b: float, 
    poly_coeffs: List[float], 
    n_points: int
) -> float:
    """
    Estime l'intgrale de P(x) sur [a,b] avec la mthode de Monte Carlo.

    Args:
        a, b: Bornes de l'intervalle.
        poly_coeffs: Coefficients du polynme, de a_0   a_n.
        n_points: Nombre de points pour le tirage.

    Returns:
        Valeur estime de l'intgrale.
    """
    if n_points <= 0:
        raise ValueError("Le nombre de points (n_points) doit tre positif.")

    # 1. Dterminer le rectangle englobant
    
    # a. Trouver les points o la drive s'annule
    derivative_coeffs = calculate_derivative_coeffs(poly_coeffs)
    derivative_roots = find_real_roots(derivative_coeffs)

    # b. Identifier les points critiques pour les extrema sur [a,b]
    critical_points = {a, b}
    for root in derivative_roots:
        if a <= root <= b:
            critical_points.add(root)

    # c. valuer le polynme aux points critiques
    extreme_values = [evaluate_polynomial(poly_coeffs, p) for p in critical_points]

    # d. Dterminer les bornes verticales du rectangle
    y_min_f = min(extreme_values)
    y_max_f = max(extreme_values)
    y_min_rect = min(0, y_min_f)
    y_max_rect = max(0, y_max_f)

    rect_area = (b - a) * (y_max_rect - y_min_rect)
    if rect_area == 0:
        return 0.0

    # 2. Tirage de Monte Carlo
    
    contributions_sum = 0
    for _ in range(n_points):
        # Tirage alatoire uniforme dans le rectangle
        rand_x = random.uniform(a, b)
        rand_y = random.uniform(y_min_rect, y_max_rect)

        f_value = evaluate_polynomial(poly_coeffs, rand_x)

        # Calcul de la contribution du point
        if 0 <= rand_y <= f_value:
            contributions_sum += 1
        elif f_value <= rand_y < 0:
            contributions_sum -= 1
            
    # 3. Calcul final
    
    estimated_integral = (contributions_sum / n_points) * rect_area
    return estimated_integral

# --- Section d'exemple ---

if __name__ == '__main__':
    # Exemple d'utilisation avec le polynme f(x) = x^2 - 4
    # Intgrale de 0   3.
    # L'intgrale exacte est [x^3/3 - 4x] de 0   3 = (27/3 - 12) - 0 = 9 - 12 = -3.
    
    # Coefficients pour x^2 - 4: a0=-4, a1=0, a2=1
    coeffs = [-4.0, 0.0, 1.0] 
    interval_a = 0.0
    interval_b = 3.0
    num_points = 100000  # Plus le nombre est grand, plus la prcision augmente

    print(f"Polynme: P(x) = x^2 - 4")
    print(f"Intervalle d'intgration: [{interval_a}, {interval_b}]")
    print(f"Nombre de points pour l'estimation: {num_points}")
    
    # Appel de la fonction principale
    estimated_value = estimate_integral_monte_carlo(
        a=interval_a,
        b=interval_b,
        poly_coeffs=coeffs,
        n_points=num_points
    )
    
    print(f"\nValeur exacte de l'intgrale: -3.0")
    print(f"Valeur estime par Monte Carlo: {estimated_value:.4f}")

    # Deuxime exemple: f(x) = -x + 5 sur [0, 5]
    # L'intgrale exacte est [-x^2/2 + 5x] de 0   5 = (-25/2 + 25) = 12.5
    coeffs_2 = [5.0, -1.0]
    interval_a_2 = 0.0
    interval_b_2 = 5.0
    
    print("\n" + "="*30)
    print(f"\nPolynme: P(x) = -x + 5")
    print(f"Intervalle d'intgration: [{interval_a_2}, {interval_b_2}]")
    
    estimated_value_2 = estimate_integral_monte_carlo(
        a=interval_a_2,
        b=interval_b_2,
        poly_coeffs=coeffs_2,
        n_points=num_points
    )
    
    print(f"\nValeur exacte de l'intgrale: 12.5")
    print(f"Valeur estime par Monte Carlo: {estimated_value_2:.4f}")
