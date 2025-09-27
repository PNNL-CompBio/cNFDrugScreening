"""
Evaluates PLSR ability to predict drug response from phosphoproteomics.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.preprocessing import normalize

from cnf.data_import import import_drug_response, import_merged_phospho
from cnf.decomposition import run_plsr
from cnf.figure_setup import get_setup


def main():
    drug_response = import_drug_response()
    merged = import_merged_phospho()

    merged = merged.loc[
        merged.index.intersection(drug_response.index),
        :
    ]
    merged.loc[:] = np.log(merged)
    merged.loc[:] = normalize(merged, axis=1)

    drug_response = drug_response.loc[merged.index, :]
    drug_response = drug_response.dropna(axis=1)

    ranks = np.arange(1, 11)
    r2s = pd.Series(0, dtype=float, index=ranks)
    fig, ax = get_setup(1, 1, {"figsize": (8, 4)})
    for rank in ranks:
        r2, predicted, plsr = run_plsr(
            merged,
            drug_response,
            n_components=int(rank)
        )
        r2s.loc[rank] = r2

    ax.plot(
        r2s.index,
        r2s,
    )
    ax.set_ylabel("R2 Score")
    ax.set_xlabel("PLSR Components")

    plt.show()


if __name__ == "__main__":
    main()
