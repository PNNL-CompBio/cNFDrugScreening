"""
Evaluates for batch effects via PCA scores plot.
"""
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize, scale

from cnf.data_import import import_merged_phospho, import_metadata


def main():
    meta = import_metadata()
    merged = import_merged_phospho()
    merged.loc[:] = normalize(merged, axis=1)
    merged.loc[:] = scale(merged, axis=0)

    cohorts = pd.Series(
        meta.loc[:, "cohort"].values,
        index=meta.loc[:, "improve_sample_id"].values
    )
    cohorts = cohorts.loc[merged.index]

    pca = PCA(n_components=2)
    components = pca.fit_transform(merged)

    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    ax.scatter(
        components[cohorts == 1, 0],
        components[cohorts == 1, 1],
        c="tab:blue",
        label="Cohort 1"
    )
    ax.scatter(
        components[cohorts == 2, 0],
        components[cohorts == 2, 1],
        c="tab:orange",
        label="Cohort 2"
    )
    ax.legend()

    plt.show()


if __name__ == "__main__":
    main()
