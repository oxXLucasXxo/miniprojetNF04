# -*- coding: utf-8 -*-
"""
Implémentation de l'intégration numérique par la méthode de Monte Carlo
pour des fonctions polynomiales.
"""

import random
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple

# --- Fonctions de l'algorithme ---

def evaluate_polynomial(coefficients: List[float], x: float) -> float:
    """
    Calcule la valeur du polynôme P(x) pour une valeur x donnée.

    Args:
        coefficients: Liste des coefficients du polynôme, de a_0 à a_n.
        x: Point d’évaluation.

    Returns:
        La valeur P(x).
    """
    # np.polyval attend les coefficients du plus haut degré au plus bas.
    return np.polyval(coefficients[::-1], x)

def calculate_derivative_coeffs(coefficients: List[float]) -> List[float]:
    """
    Calcule les coefficients du polynôme dérivé P'(x).

    Args:
        coefficients: Liste des coefficients de P(x), de a_0 à a_n.

    Returns:
        Liste des coefficients de P'(x), de a'_0 à a'_{n-1}.
    """
    if len(coefficients) < 2:
        return [0.0]

    poly_obj = np.poly1d(coefficients[::-1])
    derivative_poly_obj = np.polyder(poly_obj)
    
    return list(derivative_poly_obj.coeffs)[::-1]

def find_real_roots(coefficients: List[float]) -> List[float]:
    """
    Trouve les racines réelles d'un polynôme.

    Args:
        coefficients: Liste des coefficients du polynôme, de a_0 à a_n.

    Returns:
        Liste des racines réelles.
    """
    if not coefficients or all(c == 0 for c in coefficients):
        return []
        
    roots = np.roots(coefficients[::-1])
    return [r.real for r in roots if np.isreal(r)]

def estimate_integral_monte_carlo(
    a: float, 
    b: float, 
    poly_coeffs: List[float], 
    n_points: int
) -> Tuple[float, List[Tuple[float, float]], List[Tuple[float, float]], List[Tuple[float, float]], float, float]:
    """
    Estime l'intégrale de P(x) sur [a,b] avec la méthode de Monte Carlo.

    Args:
        a, b: Bornes de l'intervalle.
        poly_coeffs: Coefficients du polynôme, de a_0 à a_n.
        n_points: Nombre de points pour le tirage.

    Returns:
        Un tuple contenant:
        - La valeur estimée de l'intégrale.
        - Les points sous la courbe positive.
        - Les points sur la courbe négative.
        - Les points en dehors de la courbe.
        - L'ordonnée minimale du rectangle (y_min_rect).
        - L'ordonnée maximale du rectangle (y_max_rect).
    """
    if n_points <= 0:
        raise ValueError("Le nombre de points (n_points) doit être positif.")

    # 1. Déterminer le rectangle englobant
    derivative_coeffs = calculate_derivative_coeffs(poly_coeffs)
    derivative_roots = find_real_roots(derivative_coeffs)

    critical_points = {a, b}
    for root in derivative_roots:
        if a <= root <= b:
            critical_points.add(root)

    extreme_values = [evaluate_polynomial(poly_coeffs, p) for p in critical_points]

    y_min_f = min(extreme_values)
    y_max_f = max(extreme_values)
    y_min_rect = min(0, y_min_f)
    y_max_rect = max(0, y_max_f)

    rect_area = (b - a) * (y_max_rect - y_min_rect)
    if rect_area == 0:
        return 0.0, [], [], [], y_min_rect, y_max_rect

    # 2. Tirage de Monte Carlo
    contributions_sum = 0
    points_pos = []
    points_neg = []
    points_outside = []

    for _ in range(n_points):
        rand_x = random.uniform(a, b)
        rand_y = random.uniform(y_min_rect, y_max_rect)
        f_value = evaluate_polynomial(poly_coeffs, rand_x)

        if 0 <= rand_y <= f_value:
            contributions_sum += 1
            points_pos.append((rand_x, rand_y))
        elif f_value <= rand_y < 0:
            contributions_sum -= 1
            points_neg.append((rand_x, rand_y))
        else:
            points_outside.append((rand_x, rand_y))
            
    # 3. Calcul final
    estimated_integral = (contributions_sum / n_points) * rect_area
    return estimated_integral, points_pos, points_neg, points_outside, y_min_rect, y_max_rect

# --- Fonctions pour l'interaction utilisateur ---

def get_polynomial_from_user() -> List[float]:
    """Demande à l'utilisateur de saisir les coefficients du polynôme."""
    while True:
        try:
            degree_str = input("Entrez le degré du polynôme (entier >= 0): ")
            degree = int(degree_str)
            if degree < 0:
                raise ValueError("Le degré ne peut pas être négatif.")
            break
        except ValueError as e:
            print(f"Erreur: l'entrée doit être un entier positif. {e}")

    coeffs = []
    for i in range(degree + 1):
        while True:
            try:
                c_str = input(f"  - Entrez le coefficient a_{i} (pour x^{i}): ")
                coeffs.append(float(c_str))
                break
            except ValueError:
                print("Erreur: le coefficient doit être un nombre.")
    return coeffs

def get_interval_from_user() -> Tuple[float, float]:
    """Demande à l'utilisateur de saisir l'intervalle d'intégration."""
    while True:
        try:
            a_str = input("Entrez la borne inférieure de l'intervalle (a): ")
            a = float(a_str)
            b_str = input("Entrez la borne supérieure de l'intervalle (b): ")
            b = float(b_str)
            if b <= a:
                raise ValueError("La borne supérieure (b) doit être > à la borne inférieure (a).")
            return a, b
        except ValueError as e:
            print(f"Erreur: les bornes doivent être des nombres. {e}")

def get_n_points_from_user() -> int:
    """Demande à l'utilisateur le nombre de points pour l'estimation."""
    while True:
        try:
            n_str = input("Entrez le nombre de points pour l'estimation (entier > 0): ")
            n_points = int(n_str)
            if n_points <= 0:
                raise ValueError("Le nombre de points doit être un entier strictement positif.")
            return n_points
        except ValueError as e:
            print(f"Erreur: l'entrée doit être un entier positif. {e}")

def build_poly_string(coeffs: List[float]) -> str:
    """Construit une représentation textuelle du polynôme."""
    if not coeffs:
        return "0"
    parts = []
    for i, c in reversed(list(enumerate(coeffs))):
        if abs(c) < 1e-9:  # Ignorer les coefficients nuls
            continue
        
        # Signe
        sign = "-" if c < 0 else "+"
        c = abs(c)

        # Coefficient
        coeff_str = f"{c:.2f}" # Changed from {c:2f} to {c:.2f}
        
        # Terme en x
        if i == 0:
            term = coeff_str
        elif i == 1:
            term = f"{coeff_str}x" if c != 1.0 else "x"
        else:
            term = f"{coeff_str}x^{i}" if c != 1.0 else f"x^{i}"
        
        parts.append(f" {sign} {term}")

    poly_str = "".join(parts).lstrip(" +")
    return poly_str if poly_str else "0"

def plot_integral_and_points(
    a: float, 
    b: float, 
    poly_coeffs: List[float], 
    estimated_value: float,
    points_pos: List[Tuple[float, float]],
    points_neg: List[Tuple[float, float]],
    points_outside: List[Tuple[float, float]],
    y_min_rect: float,
    y_max_rect: float
):
    """
    Affiche le graphique de la fonction, l'aire d'intégration et les points de Monte Carlo.

    Args:
        a, b: Bornes de l'intervalle.
        poly_coeffs: Coefficients du polynôme.
        estimated_value: Valeur estimée de l'intégrale.
        points_pos: Points sous la courbe positive.
        points_neg: Points sur la courbe négative.
        points_outside: Points en dehors de la courbe.
        y_min_rect: Ordonnée minimale du rectangle.
        y_max_rect: Ordonnée maximale du rectangle.
    """
    n_points = len(points_pos) + len(points_neg) + len(points_outside)
    x_curve = np.linspace(a, b, 400)
    y_curve = evaluate_polynomial(poly_coeffs, x_curve)

    plt.figure(figsize=(12, 8))

    # Points de Monte Carlo
    if n_points > 0:
        # Unpack and plot points
        if points_outside:
            x_out, y_out = zip(*points_outside)
            plt.scatter(x_out, y_out, color='gray', s=1, label=f'Points extérieurs ({len(points_outside)})')
        if points_pos:
            x_pos, y_pos = zip(*points_pos)
            plt.scatter(x_pos, y_pos, color='green', s=1, label=f'Points positifs ({len(points_pos)})')
        if points_neg:
            x_neg, y_neg = zip(*points_neg)
            plt.scatter(x_neg, y_neg, color='red', s=1, label=f'Points négatifs ({len(points_neg)})')

    # Courbe de la fonction
    plt.plot(x_curve, y_curve, 'b-', linewidth=2, label=f"P(x) = {build_poly_string(poly_coeffs)}")
    
    # Aire d'intégration
    plt.fill_between(x_curve, 0, y_curve, where=y_curve>=0, facecolor='cyan', alpha=0.5, label="Aire positive")
    plt.fill_between(x_curve, 0, y_curve, where=y_curve<=0, facecolor='magenta', alpha=0.5, label="Aire négative")
    
    # Rectangle englobant
    rect = plt.Rectangle((a, y_min_rect), b - a, y_max_rect - y_min_rect,
                         edgecolor='orange', facecolor='none', linestyle='--', linewidth=2,
                         label='Rectangle englobant')
    plt.gca().add_patch(rect)

    # Ligne de l'axe x
    plt.axhline(0, color='black', linewidth=0.8)

    # Configuration du graphique
    plt.title(f"Intégration de P(x) sur [{a}, {b}] avec {n_points} points\nValeur estimée ≈ {estimated_value:.4f}")
    plt.xlabel("x")
    plt.ylabel("P(x)")
    plt.grid(True)
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.show()

# --- Section principale ---

if __name__ == '__main__':
    print("="*50)
    print("Calcul d'intégrale par la méthode de Monte Carlo")
    print("="*50)

    try:
        poly_coeffs = get_polynomial_from_user()
        interval_a, interval_b = get_interval_from_user()
        num_points = get_n_points_from_user()

        print("\n--- Récapitulatif ---")
        print(f"Polynôme P(x) = {build_poly_string(poly_coeffs)}")
        print(f"Intervalle d'intégration: [{interval_a}, {interval_b}]")
        print(f"Nombre de points pour l'estimation: {num_points:,}")
        
        # Appel de la fonction principale
        estimated_value, points_pos, points_neg, points_outside, y_min_rect, y_max_rect = estimate_integral_monte_carlo(
            a=interval_a,
            b=interval_b,
            poly_coeffs=poly_coeffs,
            n_points=num_points
        )
        
        print("\n--- Résultat ---")
        print(f"Valeur estimée de l'intégrale: {estimated_value:.6f}")

        while True:
            show_plot = input("Voulez-vous afficher le graphique de la fonction et de l'aire d'intégration ? (o/n) : ").lower()
            if show_plot in ['o', 'n']:
                break
            else:
                print("Veuillez répondre par 'o' pour oui ou 'n' pour non.")

        if show_plot == 'o':
            plot_integral_and_points(
                a=interval_a,
                b=interval_b,
                poly_coeffs=poly_coeffs,
                estimated_value=estimated_value,
                points_pos=points_pos,
                points_neg=points_neg,
                points_outside=points_outside,
                y_min_rect=y_min_rect,
                y_max_rect=y_max_rect
            )

    except ValueError as e:
        print(f"\nErreur de saisie: {e}")
    except KeyboardInterrupt:
        print("\n\nProgramme interrompu par l'utilisateur.")
    except Exception as e:
        print(f"\nUne erreur inattendue est survenue: {e}")