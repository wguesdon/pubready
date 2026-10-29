#!/usr/bin/env python3
# Standalone reproduction of the transition counts and the test. Run inside the
# pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from statsmodels.stats.contingency_tables import SquareTable

raw = pd.read_csv("input_response.csv")
d = raw[["patient", "week", "response"]].dropna().astype(str)
d.columns = ["id", "visit", "state"]
d["visit"] = d["visit"].str.replace(r"\.0$", "", regex=True)
visits = ['8', '16', '24']
fig, ax = plt.subplots(figsize=(10.09, 4.75))
for k, (a, b) in enumerate(zip(visits[:-1], visits[1:])):
    m = d[d["visit"] == a].merge(d[d["visit"] == b], on="id")
    tab = pd.crosstab(m["state_x"], m["state_y"])
    print(tab)
    if len(visits) == 2 and list(tab.index) == list(tab.columns):
        print(SquareTable(tab.to_numpy(), shift_zeros=False).symmetry(method="bowker"))
    tab.sum(axis=1).plot.bar(ax=ax, position=k, width=0.3)
ax.set_xlabel("response"); ax.set_ylabel("Patients")
fig.tight_layout(); fig.savefig("figure_response_sankey_py.pdf")
print("Reproduced the transition counts.")
