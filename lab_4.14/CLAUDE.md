# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a physics lab report project for Lab 4.14 "Spectrometer" (Лабораторная работа 4.14 «Спектрометр»). The project processes experimental measurements to calculate physical quantities: wavelengths from a diffraction grating and refractive indices from a prism.

**Main deliverable**: A LaTeX report (`report_4_14.tex`) compiled with generated plots and calculated tables.

## Project Structure

- `report_4_14.tex` — LaTeX report file (compiles with figures from `figures/`)
- `scripts/calculate_and_plot_4_14.py` — Main Python script that processes all experimental data
- `scripts/requirements.txt` — Python dependencies (numpy, pandas, matplotlib)
- `data/` — Input CSV files with raw experimental measurements
  - `grating_measurements.csv` — Diffraction grating measurements
  - `prism_measurements.csv` — Prism measurements
- `output/` — Generated results (CSV tables with calculations)
- `figures/` — Generated plots and diagrams (must be in same directory as `report_4_14.tex` for LaTeX to find them)

## Common Commands

### Setup and Run Calculations

```bash
# Create virtual environment and install dependencies
python -m venv .venv
source .venv/bin/activate  # macOS/Linux
pip install -r scripts/requirements.txt

# Run the main calculation script
python scripts/calculate_and_plot_4_14.py
```

This generates:
- `output/grating_results.csv` — Calculated wavelengths with error analysis
- `output/prism_results.csv` — Calculated refractive indices with error analysis
- `figures/` — Updated plots for LaTeX

### LaTeX Report Compilation

Compile `report_4_14.tex` in Overleaf or locally with pdflatex/xelatex. Ensure `figures/` directory is in the same location as the `.tex` file.

## Code Architecture

### Main Script (`calculate_and_plot_4_14.py`)

The script follows a straightforward structure:

1. **Configuration section** — Experiment parameters:
   - Diffraction grating: 1200 lines/mm, first order (k=1)
   - Equilateral prism apex angle: 60°
   - Ruler measurement uncertainty: ±0.1 cm for both x and y

2. **Calculation functions**:
   - `calc_grating()` — For each measurement, calculates wavelength λ from diffraction angle using `sin(θ) = λ / (k × grating spacing)`, includes error propagation
   - `calc_prism()` — For each measurement, calculates refractive index n using prism formula with minimum deviation condition, includes error propagation

3. **Data flow**: Load CSV → apply row-wise calculations → append results → save to output CSV and generate plots

### Key Physics Formulas in Code

- **Diffraction grating**: `λ_nm = 1e6 × y / (k × g × r)` where r = √(x² + y²)
- **Prism refractive index**: `n = sin((A + γ)/2) / sin(A/2)` for equilateral prism at minimum deviation

## Important Notes on Data

**Prism measurements limitation**: The original prism measurements lack color annotations for spectral lines. Lines are conditionally matched to two characteristic mercury lines (436 nm and 546 nm). If precise color information becomes available, update `data/prism_measurements.csv` and rerun the script.

## Development Notes

- All file paths use `pathlib.Path` (cross-platform compatible)
- Error propagation calculations use partial derivatives for uncertainty
- Matplotlib generates PNG figures suitable for LaTeX inclusion
- CSV output includes both calculated values and their uncertainties