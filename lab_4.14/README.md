# Лабораторная работа 4.14 «Спектрометр»

В папке лежат материалы для отчета:

- `report_4_14.tex` — готовый LaTeX-отчет.
- `data/grating_measurements.csv` — исходные измерения для дифракционной решетки.
- `data/prism_measurements.csv` — исходные измерения для призмы.
- `scripts/calculate_and_plot_4_14.py` — расчет всех таблиц и построение графиков.
- `scripts/requirements.txt` — зависимости Python.
- `output/grating_results.csv` — рассчитанная таблица для решетки.
- `output/prism_results.csv` — рассчитанная таблица для призмы.
- `figures/` — графики и схема для вставки в LaTeX.

## Как запустить расчеты локально

Из корня папки `lab_4_14_package`:

```bash
python -m venv .venv
source .venv/bin/activate  # macOS / Linux
pip install -r scripts/requirements.txt
python scripts/calculate_and_plot_4_14.py
```

После запуска обновятся файлы в папках `output/` и `figures/`.

## Как собрать отчет

Открой `report_4_14.tex` в Overleaf или локально скомпилируй через LaTeX.
Важно: папка `figures/` должна лежать рядом с `report_4_14.tex`, иначе LaTeX не найдет картинки.

## Важное замечание по данным

В исходных измерениях для призмы не были отдельно записаны цвета спектральных линий. Поэтому в отчете строки условно сопоставлены с двумя характерными линиями ртути: 436 нм и 546 нм. Если появятся точные цвета линий, достаточно заменить значения в `data/prism_measurements.csv` и заново запустить скрипт.
