"""
Runs PCA on single-cell image features.
"""
from cnf.data_import import import_synapse_folder
from cnf.figure_setup import get_setup

import pandas as pd
from sklearn.preprocessing import normalize, scale
from sklearn.decomposition import PCA
from sklearn.impute import KNNImputer
import matplotlib.pyplot as plt


def main():
    sc_data = import_synapse_folder().reset_index(drop=True)

    meta = sc_data.loc[:, :"Area.Size.Shape_Cell_CENTER.X"].iloc[:, :-1]
    sc_data = sc_data.loc[:, "Area.Size.Shape_Cell_CENTER.X":]
    sc_data = sc_data.dropna(axis=1, how="all")
    sc_data = sc_data.dropna(axis=0, how="all")
    meta = meta.loc[sc_data.index, :]
    meta.loc[:, "Metadata_Sarcoma"] = meta.loc[
        :,
        "Metadata_patient"
    ].str.contains("SARCO")

    knn = KNNImputer()
    sc_data[:] = knn.fit_transform(sc_data)
    sc_data[:] = normalize(sc_data, axis=1)
    sc_data[:] = scale(sc_data, axis=0)
    pca = PCA(n_components=10)
    factors = pd.DataFrame(
        pca.fit_transform(sc_data),
        index=sc_data.index,
        columns=range(1, pca.n_components + 1)
    )

    meta_cols = [
        "Metadata_dose", "Metadata_Well", "Metadata_patient",
        "Metadata_Target", "Metadata_Therapeutic_Categories",
        "Metadata_Sarcoma"
    ]

    fig, axes = get_setup(2, 3, {"figsize": (2 * len(meta_cols), 4)})
    axes = axes.flatten()
    for meta_col, ax in zip(meta_cols, axes):
        for group in meta.loc[:, meta_col].unique():
            ax.scatter(
                factors.loc[meta.loc[:, meta_col] == group, 1],
                factors.loc[meta.loc[:, meta_col] == group, 2],
                label=group,
                s=5
            )

        ax.set_title(meta_col)
        ax.set_xlabel(
            "PCA 1"
            f"({round(pca.explained_variance_ratio_[0] * 100, 2)}%)"
        )
        ax.set_ylabel(
            "PCA 2"
            f"({round(pca.explained_variance_ratio_[1] * 100, 2)}%)"
        )

    plt.show()


if __name__ == "__main__":
    main()
