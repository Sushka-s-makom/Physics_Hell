"""
Лабораторная работа 4.05: интерферометр Майкельсона.
Скрипт выполняет обработку измерений, строит графики и сохраняет таблицы результатов.

Запуск из корня папки проекта:
    python scripts/processing.py
"""
from __future__ import annotations

from pathlib import Path
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
FIG_DIR = ROOT / "figures_4.05"
RES_DIR = ROOT / "results"
for d in [FIG_DIR, RES_DIR]:
    d.mkdir(parents=True, exist_ok=True)

C = 2.99792458e8  # м/с
T_95_DF3 = 3.182  # для 4 значений длины волны после исключения нулевой строки
T_95_DF5 = 2.571  # для 6 значений длины когерентности
LAMBDA_TAB_NM = 546.0

# -----------------------------------------------------------------------------
# 1. Исходные данные
# -----------------------------------------------------------------------------
# В первой строке d1=d2=4, n=17: это начальное состояние, а не интервал
# перемещения. Для расчета lambda по формуле lambda = 2D/n эта строка не
# используется, так как D=0.
wavelength_measurements = pd.DataFrame(
    {
        "N": [1, 2, 3, 4, 5],
        "d1_um": [4.00, 4.00, 4.00, 4.00, 4.00],
        "d2_um": [4.00, 4.05, 4.10, 4.15, 4.20],
        "n": [17, 16, 16, 15, 14],
    }
)
wavelength_measurements["D_um"] = (wavelength_measurements["d2_um"] - wavelength_measurements["d1_um"]).abs()
wavelength_measurements["lambda_um"] = np.where(
    wavelength_measurements["D_um"] > 0,
    2 * wavelength_measurements["D_um"] / wavelength_measurements["n"],
    np.nan,
)
wavelength_measurements["lambda_nm"] = wavelength_measurements["lambda_um"] * 1000

lambda_valid = wavelength_measurements.dropna(subset=["lambda_nm"]).copy()
lambda_mean_nm = float(lambda_valid["lambda_nm"].mean())
lambda_std_nm = float(lambda_valid["lambda_nm"].std(ddof=1))
lambda_conf_nm = float(T_95_DF3 * lambda_std_nm / math.sqrt(len(lambda_valid)))
lambda_rel_error_to_tab = abs(lambda_mean_nm - LAMBDA_TAB_NM) / LAMBDA_TAB_NM * 100

# Данные по исчезновению интерференционной картины.
# d0 = 4 мкм - положение нулевой разности хода по записи в таблице.
d0_um = 4.00
coh_measurements = pd.DataFrame(
    {
        "direction": ["в большую сторону"] * 3 + ["в меньшую сторону"] * 3,
        "N": [1, 2, 3, 1, 2, 3],
        "d0_um": [d0_um] * 6,
        "d_disappear_um": [5.8, 5.9, 6.0, 3.15, 1.5, 1.7],
    }
)
coh_measurements["D_um"] = (coh_measurements["d_disappear_um"] - coh_measurements["d0_um"]).abs()
# В интерферометре Майкельсона оптическая разность хода равна удвоенному
# перемещению зеркала. Поэтому L_coh = 2D.
coh_measurements["L_coh_um"] = 2 * coh_measurements["D_um"]

L_coh_mean_um = float(coh_measurements["L_coh_um"].mean())
L_coh_std_um = float(coh_measurements["L_coh_um"].std(ddof=1))
L_coh_conf_um = float(T_95_DF5 * L_coh_std_um / math.sqrt(len(coh_measurements)))

tau_coh_s = L_coh_mean_um * 1e-6 / C
tau_coh_fs = tau_coh_s * 1e15
# Для физически корректной оценки ширины линии используем табличное значение
# зеленой линии ртути, так как экспериментальная серия для lambda записана с
# существенной ошибкой подсчета/фиксации полос.
delta_lambda_nm = LAMBDA_TAB_NM**2 / (L_coh_mean_um * 1000)  # lambda^2 / L, все в нм
# propagation: Delta(lambda_width)/lambda_width = Delta(L)/L
# (погрешностью табличной lambda пренебрегаем)
delta_lambda_unc_nm = delta_lambda_nm * L_coh_conf_um / L_coh_mean_um

delta_nu_Hz = C * (delta_lambda_nm * 1e-9) / (LAMBDA_TAB_NM * 1e-9) ** 2
delta_nu_THz = delta_nu_Hz / 1e12
delta_nu_unc_THz = delta_nu_THz * delta_lambda_unc_nm / delta_lambda_nm

# Лазер: в таблице число колец не менялось.
laser_measurements = pd.DataFrame(
    {
        "N": [1, 2, 3],
        "d1_um": [4.00, 4.00, 4.00],
        "d2_um": [4.00, 4.05, 5.60],
        "n_observed": [13, 13, 13],
    }
)
laser_measurements["D_um"] = (laser_measurements["d2_um"] - laser_measurements["d1_um"]).abs()
laser_measurements["optical_path_change_um"] = 2 * laser_measurements["D_um"]
laser_min_coh_um = float(laser_measurements["optical_path_change_um"].max())

# -----------------------------------------------------------------------------
# 2. Сохранение таблиц
# -----------------------------------------------------------------------------
wavelength_measurements.to_csv(RES_DIR / "wavelength_measurements.csv", index=False)
coh_measurements.to_csv(RES_DIR / "coherence_measurements.csv", index=False)
laser_measurements.to_csv(RES_DIR / "laser_measurements.csv", index=False)

summary = {
    "lambda_mean_raw_nm": lambda_mean_nm,
    "lambda_conf_raw_nm": lambda_conf_nm,
    "lambda_rel_error_to_tab_percent": lambda_rel_error_to_tab,
    "lambda_tab_nm_used_for_bandwidth": LAMBDA_TAB_NM,
    "L_coh_mean_um": L_coh_mean_um,
    "L_coh_conf_um": L_coh_conf_um,
    "tau_coh_fs": tau_coh_fs,
    "delta_lambda_nm": delta_lambda_nm,
    "delta_lambda_unc_nm": delta_lambda_unc_nm,
    "delta_nu_THz": delta_nu_THz,
    "delta_nu_unc_THz": delta_nu_unc_THz,
    "laser_min_coherence_um": laser_min_coh_um,
}
summary_df = pd.DataFrame([summary])
summary_df.to_csv(RES_DIR / "summary.csv", index=False)

with open(RES_DIR / "calculation_summary.txt", "w", encoding="utf-8") as f:
    f.write("Лабораторная работа 4.05: интерферометр Майкельсона\n")
    f.write("\n1. Расчет длины волны по записанной серии:\n")
    f.write(f"lambda_raw = ({lambda_mean_nm:.2f} ± {lambda_conf_nm:.2f}) нм\n")
    f.write(f"Отклонение от табличного значения 546 нм: {lambda_rel_error_to_tab:.1f}%\n")
    f.write("Вывод: эта серия записана некорректно для количественного определения lambda.\n")
    f.write("\n2. Когерентность ртутной лампы с зеленым фильтром:\n")
    f.write(f"L_coh = ({L_coh_mean_um:.2f} ± {L_coh_conf_um:.2f}) мкм\n")
    f.write(f"tau_coh = {tau_coh_fs:.2f} фс\n")
    f.write(f"Delta lambda = ({delta_lambda_nm:.1f} ± {delta_lambda_unc_nm:.1f}) нм\n")
    f.write(f"Delta nu = ({delta_nu_THz:.1f} ± {delta_nu_unc_THz:.1f}) ТГц\n")
    f.write("\n3. Лазер:\n")
    f.write(f"В доступном диапазоне разности хода интерференционная картина не исчезала.\n")
    f.write(f"Следовательно, L_coh(laser) > {laser_min_coh_um:.2f} мкм для данной серии измерений.\n")

# -----------------------------------------------------------------------------
# 3. Графики
# -----------------------------------------------------------------------------
plt.figure(figsize=(7, 4.8))
plt.plot(lambda_valid["N"], lambda_valid["lambda_nm"], marker="o", label=r"$\lambda_N=2D/n$")
plt.axhline(LAMBDA_TAB_NM, linestyle="--", label=r"$\lambda_{tab}=546$ нм")
plt.xlabel("Номер измерения")
plt.ylabel(r"$\lambda_N$, нм")
plt.title("Расчет длины волны по записанной серии")
plt.grid(True, alpha=0.35)
plt.legend()
plt.tight_layout()
plt.savefig(FIG_DIR / "wavelength_raw.png", dpi=300)
plt.close()

plt.figure(figsize=(7, 4.8))
x_labels = [f"{d}\n{n}" for d, n in zip(coh_measurements["direction"], coh_measurements["N"])]
plt.bar(x_labels, coh_measurements["L_coh_um"])
plt.axhline(L_coh_mean_um, linestyle="--", label=rf"среднее {L_coh_mean_um:.2f} мкм")
plt.ylabel(r"$L_{ког}$, мкм")
plt.title("Оценка длины когерентности по исчезновению колец")
plt.xticks(rotation=20, ha="right")
plt.grid(axis="y", alpha=0.35)
plt.legend()
plt.tight_layout()
plt.savefig(FIG_DIR / "coherence_length.png", dpi=300)
plt.close()

plt.figure(figsize=(7, 4.8))
plt.plot(laser_measurements["d2_um"], laser_measurements["n_observed"], marker="o")
plt.xlabel(r"Положение винта $d_2$, мкм")
plt.ylabel("Число наблюдаемых колец")
plt.title("Наблюдение картины при использовании лазера")
plt.grid(True, alpha=0.35)
plt.tight_layout()
plt.savefig(FIG_DIR / "laser_observation.png", dpi=300)
plt.close()

print("Готово. Результаты сохранены в:")
print(f"  {RES_DIR}")
print(f"  {FIG_DIR}")
