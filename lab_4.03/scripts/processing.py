from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "data" / "measurements_clean.csv"
FIG_DIR = ROOT / "figure_4.03"
RES_DIR = ROOT / "results"
FIG_DIR.mkdir(exist_ok=True)
RES_DIR.mkdir(exist_ok=True)

# Константы обработки
T_95_DF7 = 2.365  # коэффициент Стьюдента для P=0.95, df=7
DR_DISAPPEAR_UM = 25.0  # погрешность визуального определения точки нулевой видности
C_M_PER_S = 2.99792458e8

FILE_NAMES = {
    "Зелёный": "green",
    "Красный": "red",
    "Оранжевый": "orange",
    "Синий": "blue",
}


def nm_to_cm(value_nm: float) -> float:
    return value_nm * 1e-7


def um_to_cm(value_um: float) -> float:
    return value_um * 1e-4


def calc_bandwidth(lambda_nm: float, r_disappear_um: float, R_cm: float) -> float:
    """Ширина пропускания фильтра по формуле из методички.

    Все вычисления внутри функции выполняются в сантиметрах.
    Возвращается результат в нанометрах.
    """
    lam_cm = nm_to_cm(lambda_nm)
    r_cm = um_to_cm(r_disappear_um)
    delta_lambda_cm = 2 * lam_cm**2 * R_cm / (2 * r_cm**2 + R_cm * lam_cm)
    return delta_lambda_cm / 1e-7


def calc_bandwidth_uncertainty(lambda_nm: float, r_disappear_um: float, R_cm: float,
                              dR_cm: float, dr_um: float = DR_DISAPPEAR_UM) -> float:
    """Погрешность ширины пропускания через частные производные.

    Считаем lambda заданной точно, а погрешности вносят R и r_disappear.
    Возвращается результат в нанометрах.
    """
    lam = nm_to_cm(lambda_nm)
    r = um_to_cm(r_disappear_um)
    dr = um_to_cm(dr_um)
    denom = 2 * r**2 + R_cm * lam
    # f = 2 lambda^2 R / (2 r^2 + R lambda)
    df_dR = (4 * lam**2 * r**2) / denom**2
    df_dr = -(8 * lam**2 * R_cm * r) / denom**2
    d_delta_lambda_cm = math.sqrt((df_dR * dR_cm) ** 2 + (df_dr * dr) ** 2)
    return d_delta_lambda_cm / 1e-7


def main() -> None:
    df = pd.read_csv(DATA_PATH)
    radii_rows = []
    radius_rows = []
    all_R_values = []

    for filter_name, g in df.groupby("filter", sort=False):
        wavelength_nm = float(g["wavelength_nm"].iloc[0])
        wavelength_cm = nm_to_cm(wavelength_nm)
        rings = g["ring"].to_numpy(dtype=float)

        diameters_um = g[["diameter_1_um", "diameter_2_um", "diameter_3_um"]].to_numpy(dtype=float)
        radii_um = diameters_um / 2.0
        r_mean_um = radii_um.mean(axis=1)
        r_std_um = radii_um.std(axis=1, ddof=1)
        r_mean_cm = um_to_cm(r_mean_um)
        r2_cm2 = r_mean_cm**2

        # МНК: r^2 = k n + b
        A = np.vstack([rings, np.ones_like(rings)]).T
        slope_cm2, intercept_cm2 = np.linalg.lstsq(A, r2_cm2, rcond=None)[0]
        R_mnk_cm = slope_cm2 / wavelength_cm

        # Расчёт по двум парам колец: (4,2) и (3,1)
        def R_pair(m_ring: int, n_ring: int) -> float:
            m_idx = np.where(rings == m_ring)[0][0]
            n_idx = np.where(rings == n_ring)[0][0]
            return (r_mean_cm[m_idx] ** 2 - r_mean_cm[n_idx] ** 2) / ((m_ring - n_ring) * wavelength_cm)

        R_42_cm = R_pair(4, 2)
        R_31_cm = R_pair(3, 1)
        R_avg_filter_cm = (R_42_cm + R_31_cm) / 2
        all_R_values.extend([R_42_cm, R_31_cm])

        for i, ring in enumerate(rings.astype(int)):
            radii_rows.append({
                "filter": filter_name,
                "wavelength_nm": wavelength_nm,
                "ring": ring,
                "r1_um": radii_um[i, 0],
                "r2_um": radii_um[i, 1],
                "r3_um": radii_um[i, 2],
                "r_mean_um": r_mean_um[i],
                "r_std_um": r_std_um[i],
                "r2_cm2": r2_cm2[i],
            })

        radius_rows.append({
            "filter": filter_name,
            "wavelength_nm": wavelength_nm,
            "slope_cm2_per_order": slope_cm2,
            "intercept_cm2": intercept_cm2,
            "R_MNK_cm": R_mnk_cm,
            "R_42_cm": R_42_cm,
            "R_31_cm": R_31_cm,
            "R_avg_filter_cm": R_avg_filter_cm,
        })

        # График r^2(n) для каждого фильтра
        x_line = np.linspace(1, 4, 100)
        y_line = slope_cm2 * x_line + intercept_cm2
        plt.figure(figsize=(6.5, 4.4))
        plt.scatter(rings, r2_cm2, label="Эксперимент")
        plt.plot(x_line, y_line, label="МНК-аппроксимация")
        plt.xlabel("Номер тёмного кольца n")
        plt.ylabel(r"$r^2$, см$^2$")
        plt.title(f"{filter_name} светофильтр, λ = {wavelength_nm:g} нм")
        plt.grid(True, alpha=0.35)
        plt.legend()
        plt.tight_layout()
        plt.savefig(FIG_DIR / f"r2_{FILE_NAMES[filter_name]}.png", dpi=300)
        plt.close()

    all_R_values = np.array(all_R_values, dtype=float)
    R_mean_cm = float(all_R_values.mean())
    R_std_cm = float(all_R_values.std(ddof=1))
    R_conf_cm = float(T_95_DF7 * R_std_cm / math.sqrt(len(all_R_values)))

    # Полоса пропускания фильтров
    bandwidth_rows = []
    for filter_name, g in df.groupby("filter", sort=False):
        wavelength_nm = float(g["wavelength_nm"].iloc[0])
        zero_values = g[["zero_visibility_1_um", "zero_visibility_2_um", "zero_visibility_3_um"]].iloc[0].to_numpy(dtype=float)
        r_disappear_mean_um = float(np.mean(zero_values))
        r_disappear_std_um = float(np.std(zero_values, ddof=1))
        delta_lambda_nm = calc_bandwidth(wavelength_nm, r_disappear_mean_um, R_mean_cm)
        delta_lambda_unc_nm = calc_bandwidth_uncertainty(wavelength_nm, r_disappear_mean_um, R_mean_cm, R_conf_cm)
        # Частотная полоса: |Delta nu| = c * Delta lambda / lambda^2
        delta_nu_THz = C_M_PER_S * (delta_lambda_nm * 1e-9) / (wavelength_nm * 1e-9) ** 2 / 1e12
        delta_nu_unc_THz = C_M_PER_S * (delta_lambda_unc_nm * 1e-9) / (wavelength_nm * 1e-9) ** 2 / 1e12

        bandwidth_rows.append({
            "filter": filter_name,
            "wavelength_nm": wavelength_nm,
            "r_disappear_1_um": zero_values[0],
            "r_disappear_2_um": zero_values[1],
            "r_disappear_3_um": zero_values[2],
            "r_disappear_mean_um": r_disappear_mean_um,
            "r_disappear_std_um": r_disappear_std_um,
            "adopted_dr_disappear_um": DR_DISAPPEAR_UM,
            "delta_lambda_nm": delta_lambda_nm,
            "delta_lambda_unc_nm": delta_lambda_unc_nm,
            "delta_nu_THz": delta_nu_THz,
            "delta_nu_unc_THz": delta_nu_unc_THz,
        })

    pd.DataFrame(radii_rows).to_csv(RES_DIR / "radii_summary.csv", index=False)
    pd.DataFrame(radius_rows).to_csv(RES_DIR / "radius_results.csv", index=False)
    pd.DataFrame(bandwidth_rows).to_csv(RES_DIR / "bandwidth_results.csv", index=False)

    with open(RES_DIR / "calculation_summary.txt", "w", encoding="utf-8") as f:
        f.write("Лабораторная работа 4.03. Кольца Ньютона\n")
        f.write("\nИтоговый радиус кривизны линзы:\n")
        f.write(f"R = ({R_mean_cm:.1f} ± {R_conf_cm:.1f}) см\n")
        f.write(f"R = ({R_mean_cm/100:.3f} ± {R_conf_cm/100:.3f}) м\n")
        f.write("\nРезультаты по светофильтрам:\n")
        for row in bandwidth_rows:
            f.write(
                f"{row['filter']}: lambda = {row['wavelength_nm']} нм, "
                f"r_disappear = {row['r_disappear_mean_um']:.4f} мкм, "
                f"Delta lambda = ({row['delta_lambda_nm']:.1f} ± {row['delta_lambda_unc_nm']:.1f}) нм, "
                f"Delta nu = ({row['delta_nu_THz']:.1f} ± {row['delta_nu_unc_THz']:.1f}) ТГц\n"
            )

    print((RES_DIR / "calculation_summary.txt").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
