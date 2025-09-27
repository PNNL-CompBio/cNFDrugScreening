"""
Utilities for decomposing datasets.
"""
from typing import Tuple

import pandas as pd
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import r2_score
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler


def run_plsr(
    data: pd.DataFrame,
    target: pd.DataFrame,
    n_components: int = 2
) -> Tuple[float, pd.DataFrame, PLSRegression]:
    """
    Runs and evaluates PLSR performance.

    Args:
        data (pd.DataFrame): Data to regress.
        target (pd.DataFrame): Target data to regress against.
        n_components (int): Number of PLSR components to use.

    Returns:
        float: R2 score for regression.
        pd.DataFrame: Predicted target values via PLSR.
        PLSRegression: PLSR model.
    """
    plsr = PLSRegression(n_components=n_components, scale=False)
    loo = LeaveOneOut()
    data_scaler = StandardScaler()
    target_scaler = StandardScaler()
    predicted = pd.DataFrame(
        0,
        dtype=float,
        index=target.index,
        columns=target.columns
    )

    # for train_index, test_index in kf.split(data, target):
    for train_index, test_index in loo.split(data, target):
        train_data = data.iloc[train_index, :]
        train_target = target.iloc[train_index, :]
        test_data = data.iloc[test_index, :]

        train_data = data_scaler.fit_transform(train_data)
        train_target = target_scaler.fit_transform(train_target)
        test_data = data_scaler.transform(test_data)

        plsr.fit(train_data, train_target)
        predicted.iloc[test_index, :] = target_scaler.inverse_transform(
            plsr.predict(test_data)
        )

    r2 = r2_score(target, predicted)
    plsr.fit(
        data_scaler.fit_transform(data),
        target_scaler.fit_transform(target)
    )

    return r2, predicted, plsr
