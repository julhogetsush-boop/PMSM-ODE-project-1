# Clean reproduction note

Run from the repository root:

```powershell
python code/run_all.py
```

The run regenerates `results/summary.json`, 24 CSV files and 14 figures.
A suitable individual-challenge line is:

```text
Results from your run: Explicit Euler observed order = 1.016380216709898; standard error = 0.0034410504662547425. This rounds to 1.016380.
```

The fitted range is the four finest declared meshes: 50, 25, 12.5 and 6.25 microseconds.