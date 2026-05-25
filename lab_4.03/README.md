# Лабораторная работа 4.03

Тема: **Определение радиуса кривизны линзы по интерференционной картине колец Ньютона**.

## Структура папки

```text
lab_4_03_newton_rings/
├── main.tex                         # LaTeX-отчёт
├── report_4_03.pdf                  # собранный PDF отчёта
├── requirements.txt                 # зависимости для Python
├── data/
│   ├── measurements_original.csv     # исходный CSV
│   └── measurements_clean.csv        # аккуратно переписанные данные для обработки
├── scripts/
│   └── processing.py                 # расчёты и построение графиков
├── figure_4.03/
│   ├── r2_green.png
│   ├── r2_red.png
│   ├── r2_orange.png
│   └── r2_blue.png
└── results/
    ├── radii_summary.csv
    ├── radius_results.csv
    ├── bandwidth_results.csv
    └── calculation_summary.txt
```

## Как запустить расчёты и заново построить графики

```bash
cd lab_4.03
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/processing.py
```

## Как собрать отчёт

Отчёт рассчитан на компиляцию через XeLaTeX:

```bash
xelatex main.tex
xelatex main.tex
```

Если графики были перестроены, они автоматически сохраняются в папку `figure_4.03/`.
