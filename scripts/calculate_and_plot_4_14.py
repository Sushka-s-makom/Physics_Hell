"""
Лабораторная работа 4.14 «Спектрометр».
Скрипт рассчитывает длины волн по дифракционной решетке,
показатель преломления призмы и строит графики для отчета.

Запуск из корня папки lab_4_14_package:
    python scripts/calculate_and_plot_4_14.py

Или из папки scripts:
    python calculate_and_plot_4_14.py
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------
# Исходные настройки эксперимента
# -----------------------------
ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
OUT_DIR = ROOT / "output"
FIG_DIR = ROOT / "figures"
OUT_DIR.mkdir(exist_ok=True)
FIG_DIR.mkdir(exist_ok=True)

# Дифракционная решетка: 1200 штрихов/мм, первый порядок
G_LINES_PER_MM = 1200.0
K_ORDER = 1.0

# Погрешность линейки: 1 мм = 0.1 см
DX_CM = 0.1
DY_CM = 0.1

# Угол при вершине равносторонней призмы
A_DEG = 60.0
A_RAD = math.radians(A_DEG)


# -----------------------------
# Расчеты для дифракционной решетки
# -----------------------------
def calc_grating(row: pd.Series) -> pd.Series:
    x = float(row["x_cm"])
    y = float(row["y_cm"])
    lambda_theor = float(row["lambda_theor_nm"])

    theta_rad = math.atan2(y, x)
    theta_deg = math.degrees(theta_rad)

    # sin(theta) = y / sqrt(x^2 + y^2)
    # lambda_mm = sin(theta) / (k * g), где g в штр/мм
    # lambda_nm = lambda_mm * 1e6
    r = math.sqrt(x**2 + y**2)
    lambda_exp_nm = 1e6 * y / (K_ORDER * G_LINES_PER_MM * r)

    # Распространение погрешностей для lambda = C * y / sqrt(x^2+y^2)
    # C = 1e6 / (k*g)
    C = 1e6 / (K_ORDER * G_LINES_PER_MM)
    dlam_dx = abs(C * (-x * y) / (x**2 + y**2) ** 1.5)
    dlam_dy = abs(C * (x**2) / (x**2 + y**2) ** 1.5)
    delta_lambda_nm = math.sqrt((dlam_dx * DX_CM) ** 2 + (dlam_dy * DY_CM) ** 2)

    rel_error_percent = abs(lambda_exp_nm - lambda_theor) / lambda_theor * 100

    return pd.Series(
        {
            "theta_deg": theta_deg,
            "lambda_exp_nm": lambda_exp_nm,
            "delta_lambda_nm": delta_lambda_nm,
            "rel_error_percent": rel_error_percent,
        }
    )


# -----------------------------
# Расчеты для призмы
# -----------------------------
def calc_prism(row: pd.Series) -> pd.Series:
    x = float(row["x_cm"])
    y = float(row["y_cm"])
    n_theor = float(row["n_theor_F2"])

    gamma_rad = math.atan2(y, x)
    gamma_deg = math.degrees(gamma_rad)

    # Для равносторонней призмы в положении минимального отклонения:
    # n = sin((A + gamma)/2) / sin(A/2)
    n_exp = math.sin((A_RAD + gamma_rad) / 2) / math.sin(A_RAD / 2)

    # Погрешность gamma = arctan(y/x)
    dgamma_dx = -y / (x**2 + y**2)
    dgamma_dy = x / (x**2 + y**2)
    delta_gamma_rad = math.sqrt((dgamma_dx * DX_CM) ** 2 + (dgamma_dy * DY_CM) ** 2)
    delta_gamma_deg = math.degrees(delta_gamma_rad)

    # Погрешность n по gamma
    dn_dgamma = 0.5 * math.cos((A_RAD + gamma_rad) / 2) / math.sin(A_RAD / 2)
    delta_n = abs(dn_dgamma) * delta_gamma_rad

    # Теоретический минимальный угол отклонения для данного n_theor:
    # gamma = 2*arcsin(n*sin(A/2)) - A
    gamma_theor_rad = 2 * math.asin(n_theor * math.sin(A_RAD / 2)) - A_RAD
    gamma_theor_deg = math.degrees(gamma_theor_rad)

    rel_error_n_percent = abs(n_exp - n_theor) / n_theor * 100

    return pd.Series(
        {
            "gamma_exp_deg": gamma_deg,
            "delta_gamma_deg": delta_gamma_deg,
            "gamma_theor_deg": gamma_theor_deg,
            "n_exp": n_exp,
            "delta_n": delta_n,
            "rel_error_n_percent": rel_error_n_percent,
        }
    )


# -----------------------------
# Основной запуск
# -----------------------------
def main() -> None:
    grating = pd.read_csv(DATA_DIR / "grating_measurements.csv")
    prism = pd.read_csv(DATA_DIR / "prism_measurements.csv")

    grating_results = pd.concat([grating, grating.apply(calc_grating, axis=1)], axis=1)
    prism_results = pd.concat([prism, prism.apply(calc_prism, axis=1)], axis=1)

    # Сохраняем результаты с нормальной точностью
    grating_results.to_csv(OUT_DIR / "grating_results.csv", index=False, float_format="%.6f")
    prism_results.to_csv(OUT_DIR / "prism_results.csv", index=False, float_format="%.6f")

    # Печатаем в терминал для проверки
    pd.set_option("display.max_columns", None)
    print("\n=== Дифракционная решетка ===")
    print(grating_results.to_string(index=False))
    print("\n=== Призма ===")
    print(prism_results.to_string(index=False))

    # -----------------------------
    # График 1: сравнение lambda_exp и lambda_theor
    # -----------------------------
    fig, ax = plt.subplots(figsize=(7, 5), dpi=160)
    ax.errorbar(
        grating_results["lambda_theor_nm"],
        grating_results["lambda_exp_nm"],
        yerr=grating_results["delta_lambda_nm"],
        fmt="o",
        capsize=4,
        label="Эксперимент",
    )
    min_l = min(grating_results["lambda_theor_nm"].min(), grating_results["lambda_exp_nm"].min()) - 20
    max_l = max(grating_results["lambda_theor_nm"].max(), grating_results["lambda_exp_nm"].max()) + 20
    ax.plot([min_l, max_l], [min_l, max_l], linestyle="--", label="Идеальное совпадение")
    for _, row in grating_results.iterrows():
        ax.annotate(str(row["color"]), (row["lambda_theor_nm"], row["lambda_exp_nm"]), xytext=(6, 5), textcoords="offset points")
    ax.set_xlabel(r"$\lambda_{\mathrm{теор}}$, нм")
    ax.set_ylabel(r"$\lambda_{\mathrm{эксп}}$, нм")
    ax.set_title("Сравнение экспериментальных и табличных длин волн")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "wavelength_comparison.png")
    plt.close(fig)

    # -----------------------------
    # График 2: n(lambda) для призмы
    # -----------------------------
    fig, ax = plt.subplots(figsize=(7, 5), dpi=160)
    ax.errorbar(
        prism_results["assumed_lambda_nm"],
        prism_results["n_exp"],
        yerr=prism_results["delta_n"],
        fmt="o",
        capsize=4,
        label="Эксперимент",
    )
    ax.plot(
        prism_results["assumed_lambda_nm"],
        prism_results["n_theor_F2"],
        marker="s",
        linestyle="--",
        label="F2, табличные значения",
    )
    ax.set_xlabel(r"$\lambda$, нм")
    ax.set_ylabel(r"$n$")
    ax.set_title("Показатель преломления призмы")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "prism_refractive_index.png")
    plt.close(fig)

    # -----------------------------
    # Схематический рисунок установки с решеткой
    # -----------------------------
    fig, ax = plt.subplots(figsize=(8, 3.2), dpi=160)
    ax.set_axis_off()
    ax.set_xlim(-1, 11)
    ax.set_ylim(-2, 2)

    # Основная оптическая ось
    ax.plot([0, 10], [0, 0], linewidth=1)

    elements = [
        (0.5, "Источник\nсвета"),
        (2.2, "Линза\n$f=50$ мм"),
        (3.9, "Щель"),
        (5.7, "Линза\n$f=100$ мм"),
        (7.5, "Решетка/\nпризма"),
        (9.7, "Экран"),
    ]
    for x_pos, label in elements:
        ax.plot([x_pos, x_pos], [-0.65, 0.65], linewidth=2)
        ax.text(x_pos, -1.15, label, ha="center", va="top", fontsize=9)

    # Лучи
    ax.plot([0.5, 3.9, 7.5, 9.7], [0, 0, 0, 0.9], linestyle="-")
    ax.plot([0.5, 3.9, 7.5, 9.7], [0, 0, 0, -0.9], linestyle="-")
    ax.plot([7.5, 9.7], [0, 0], linestyle="--")

    # Подписи расстояний
    ax.annotate("24 см", xy=(2.2, 0.35), xytext=(2.2, 0.85), ha="center", arrowprops=dict(arrowstyle="<->"))
    ax.annotate("20--25 см", xy=(5.7, 0.35), xytext=(5.7, 1.25), ha="center", arrowprops=dict(arrowstyle="<->"))
    ax.annotate("$x$", xy=(8.6, -0.35), xytext=(8.6, -1.65), ha="center", arrowprops=dict(arrowstyle="<->"))
    ax.text(9.9, 0.75, "$y$", fontsize=12)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "spectrometer_scheme.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
