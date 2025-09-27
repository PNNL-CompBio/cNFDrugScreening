"""
Functions for importing data from Synapse.
"""
from os import PathLike
from os.path import abspath, dirname, join, exists

from inmoose.pycombat import pycombat_seq
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
import synapseclient as sc

REPO_PATH = dirname(dirname(abspath(__file__)))


def syn_login(auth_path: PathLike | None = None) -> sc.Synapse:
    """
    Login to Synapse.

    Args:
        auth_path (PathLike | None): Path to authentication file.

    Returns:
        sc.Synapse: Logged-in Synapse client.
    """
    if auth_path is None:
        auth_path = join(REPO_PATH, "auth_token.txt")

    syn = sc.Synapse()
    with open(auth_path, "r") as f:
        auth_token = f.read()

    syn.login(authToken=auth_token)

    return syn


def import_synapse_folder(
    syn: sc.Synapse | None = None,
    folder_id: str | None = None
) -> pd.DataFrame:
    """
    Loads all data from synapse folder.

    Args:
        syn (sc.Synapse, default: None): Logged-in Synapse client; initializes
            new one if None is provided.
        folder_id (str, default: None): Synapse folder to read; defaults to
            minimally-processed single-cell image data.

    Returns:
        pd.DataFrame: DataFrame containing concatenation of all data in folder.
    """
    if syn is None:
        syn = syn_login()

    if folder_id is None:
        folder_id = "syn69874135"

    if exists(join(REPO_PATH, "data", f"{folder_id}.parquet")):
        return pd.read_parquet(
            join(REPO_PATH, "data", f"{folder_id}.parquet")
        )

    children = syn.getChildren(parent=folder_id)
    data = None
    for data_file in children:
        if data_file["type"] == "org.sagebionetworks.repo.model.FileEntity":
            if data is None:
                data = pd.read_parquet(
                    syn.get(data_file["id"]).path
                )
                data.loc[:, "Metadata_Drug"] = data_file["name"].split("_")[2]
            else:
                _data = pd.read_parquet(
                    syn.get(data_file["id"]).path
                )
                _data.loc[:, "Metadata_Drug"] = data_file["name"].split("_")[2]
                data = pd.concat([
                    data,
                    _data
                ])

    data = data.loc[:, ["Metadata_Drug"] + list(data.columns[:-1])]
    data.to_parquet(join(REPO_PATH, "data", f"{folder_id}.parquet"))

    return data


def import_proteomics(syn: sc.Synapse | None = None) -> pd.DataFrame:
    """
    Loads proteomic data.

    Args:
        syn (sc.Synapse): Logged-in Synapse client; initializes new one if None
            is provided.

    Returns:
        pd.DataFrame: Proteomic data.
    """
    if syn is None:
        syn = syn_login()

    proteomics = pd.read_csv(
        syn.get("syn64906445").path,
        sep="\t",
        index_col=0
    )
    return proteomics


def import_phospho(
    syn: sc.Synapse | None = None
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Loads phosphoproteomic data.

    Args:
        syn (sc.Synapse): Logged-in Synapse client; initializes new one if None
            is provided.

    Returns:
        pd.DataFrame: Phosphoproteomic data for cohort 1.
    """
    if syn is None:
        syn = syn_login()

    phospho_c1 = pd.read_csv(
        syn.get("syn65467785").path,
        sep="\t",
        index_col=0
    )
    phospho_c2 = pd.read_csv(
        join("data", "DIA_Phospho_FP_Results_SiteID_c2.txt"),
        index_col=0,
        sep="\t"
    )

    return phospho_c1, phospho_c2


def import_merged_phospho():
    """
    Batch corrects and loads merged phosphoproteomic datasets.

    Args:
        use_improve (bool, default: True): Use Improve Sample IDs for sample
            names.

    Returns:
        pd.DataFrame: Merged phosphoproteomic data.
    """
    if exists(join(REPO_PATH, "data", "merged_phospho.txt.gz")):
        merged = pd.read_csv(
            join(REPO_PATH, "data", "merged_phospho.txt.gz"),
            index_col=0
        )
    else:
        phospho_1, phospho_2 = import_phospho()

        merged = pd.concat([phospho_1, phospho_2], axis=1)
        merged = merged.dropna(axis=0)
        cohorts = np.repeat(
            [1, 2],
            [phospho_1.shape[1], phospho_2.shape[1]],
            axis=0
        )
        merged = pycombat_seq(merged, cohorts)
        merged = merged.T

        meta = import_metadata()
        ids = merged.index.str.split(".", expand=True)
        ids = ["-".join(column[6:10])[8:] for column in ids]
        merged.index = ids

        sample_ids = pd.Series(
            meta.loc[:, "improve_sample_id"].values,
            index=meta.loc[:, "DatasetNamePhospho"].values
        )
        merged = merged.rename(index=sample_ids)

    # Saves batch-corrected, merged data
    merged.to_csv(join(REPO_PATH, "data", "merged_phospho.txt.gz"))

    return merged


def import_drug_response(
    syn: sc.Synapse | None = None
) -> pd.DataFrame:
    """
    Loads drug response curves.

    Args:
        syn (sc.Synapse): Logged-in Synapse client; initializes new one if None
            is provided.

    Returns:
        pd.DataFrame: Drug response curves.
    """
    if syn is None:
        syn = syn_login()

    c1 = pd.read_csv(
        syn.get("syn65471817").path,
        sep="\t",
        index_col=0
    )
    c2 = pd.read_csv(
        join(REPO_PATH, "data", "cohort2_curves.tsv"),
        sep="\t",
        index_col=0
    )

    drug_response = pd.concat([c1, c2])
    drug_response.reset_index(drop=True, inplace=True)
    drug_response = drug_response.loc[
        drug_response.loc[:, "dose_response_metric"] == "uM_viability",
        :
    ]
    drug_response = drug_response.sort_values(
        "dose_response_value",
        ascending=False
    )
    drug_response = drug_response.loc[
        ~drug_response.loc[
            :,
            ["improve_sample_id", "improve_drug_id"]
        ].duplicated(keep="first"),
        :
    ]

    drug_response = drug_response.pivot(
        index="improve_sample_id",
        columns="improve_drug_id",
        values="dose_response_value"
    )

    return drug_response


def import_metadata(syn: sc.Synapse | None = None) -> pd.DataFrame:
    """
    Loads phosphoproteomic and proteomic metadata.

    Args:
        syn (sc.Synapse): Logged-in Synapse client; initializes new one if None
            is provided.

    Returns:
        pd.DataFrame: Metadata across cohorts.
    """
    if syn is None:
        syn = syn_login()

    c1_meta = pd.read_excel(
        syn.get("syn69920463").path,
    )
    c2_meta = pd.read_excel(
        syn.get("syn69920464").path,
        sheet_name="Sheet2"
    )

    c1_meta.rename(
        columns={
            "SampleAlias": "SampleNo",
            "Count (k)": "Count"
        },
        inplace=True
    )
    c1_descriptors = c1_meta.loc[:, "Description"].str.split(" ", expand=True)
    c1_meta.loc[:, "Tumor"] = c1_descriptors.iloc[:, 1].astype(int)
    c1_meta.drop("Description.1", axis=1, inplace=True)

    c2_descriptors = c2_meta.loc[:, "Description"].str.split(
        "-",
        expand=True
    )
    c2_descriptors.columns = ["Patient", "Tumor", "Count", "Date"]
    c2_descriptors.loc[:, "Count"].replace(
        {"Proteomic 300K": 300000},
        inplace=True
    )
    c2_descriptors.loc[:, "Tumor"] = c2_descriptors.loc[
        :,
        "Tumor"
    ].str[1:].astype(int)
    c2_meta.loc[:, c2_descriptors.columns] = c2_descriptors
    c2_meta.rename(columns={"SampleAlias": "SampleNo"}, inplace=True)

    meta = pd.concat([c1_meta, c2_meta])
    meta.loc[:, "improve_sample_id"] = (
        meta.loc[:, "Patient"] +
        "_T" + meta.loc[:, "Tumor"].astype(str)
    )
    meta.loc[:, "cohort"] = np.repeat(
        [1, 2],
        [c1_meta.shape[0], c2_meta.shape[0]]
    )
    meta.loc[
        meta.loc[:, "improve_sample_id"].duplicated(),
        "improve_sample_id"
    ] += "_2"

    return meta


def reorder_table(df) -> pd.DataFrame:
    """
    Reorder a DataFrame's rows using hierarchical clustering.

    Args:
        df (pandas.DataFrame): data to be clustered; rows are treated as samples
            to be clustered.

    Returns:
        pd.DataFrame: data with rows reordered via hierarchical clustering.
    """
    y = sch.linkage(df.to_numpy(), method='centroid')
    index = sch.dendrogram(y, orientation='right', no_plot=True)['leaves']
    return df.iloc[index]


if __name__ == "__main__":
    import_phospho()
