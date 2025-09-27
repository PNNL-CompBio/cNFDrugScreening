"""
Plots heatmaps of SGCCA associations.
"""
from os.path import abspath, dirname, join

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from cnf.data_import import reorder_table
from cnf.figure_setup import get_setup

REPO_PATH = abspath(dirname(dirname(dirname(__file__))))


def main():
    p_scores = pd.read_csv(
        join(
            REPO_PATH,
            "output",
            "phospho_scores.csv"
        ),
        index_col=0
    )
    p_loadings = pd.read_csv(
        join(
            REPO_PATH,
            "output",
            "phospho_loadings.csv"
        ),
        index_col=0
    )
    dr_scores = pd.read_csv(
        join(
            REPO_PATH,
            "output",
            "drug_scores.csv"
        ),
        index_col=0
    )
    dr_loadings = pd.read_csv(
        join(
            REPO_PATH,
            "output",
            "drug_loadings.csv"
        ),
        index_col=0
    )

    fig, axes = get_setup(1, 4, {"figsize": (8, 4)})
    result_sets = [(p_scores, p_loadings), (dr_scores, dr_loadings)]
    for col_index, (scores, loadings) in enumerate(result_sets):
        scores = reorder_table(scores)
        loadings = reorder_table(loadings)
        sns.heatmap(
            scores,
            vmin=-1,
            vmax=1,
            cmap="coolwarm",
            ax=axes[col_index],
            cbar=col_index == len(result_sets) - 1
        )
        sns.heatmap(
            loadings,
            center=0,
            cmap="coolwarm",
            ax=axes[col_index + 2]
        )
        axes[col_index].set_xticklabels(np.arange(1, 5), rotation=0)
        axes[col_index + 2].set_xticklabels(np.arange(1, 5), rotation=0)

    axes[1].set_yticks([])

    plt.show()


if __name__ == "__main__":
    main()
