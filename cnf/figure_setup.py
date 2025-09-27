"""
Functions for figure setup and consistent style.
"""
from typing import Any, Dict

import matplotlib
import matplotlib.pyplot as plt
import numpy.typing as npt

matplotlib.rcParams["axes.labelsize"] = 10
matplotlib.rcParams["axes.linewidth"] = 0.6
matplotlib.rcParams["axes.titlesize"] = 12
matplotlib.rcParams["font.family"] = ["sans-serif"]
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams["font.size"] = 8
matplotlib.rcParams["grid.linestyle"] = "dotted"
matplotlib.rcParams["legend.borderpad"] = 0.35
matplotlib.rcParams["legend.fontsize"] = 7
matplotlib.rcParams["legend.framealpha"] = 0.5
matplotlib.rcParams["legend.handlelength"] = 0.5
matplotlib.rcParams["legend.handletextpad"] = 0.5
matplotlib.rcParams["legend.labelspacing"] = 0.2
matplotlib.rcParams["legend.markerscale"] = 0.7
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["xtick.labelsize"] = 8
matplotlib.rcParams["xtick.major.pad"] = 1.0
matplotlib.rcParams["xtick.minor.pad"] = 0.9
matplotlib.rcParams["ytick.labelsize"] = 8
matplotlib.rcParams["ytick.major.pad"] = 1.0
matplotlib.rcParams["ytick.minor.pad"] = 0.9


def get_setup(
    n_rows: int,
    n_cols: int,
    fig_params: Dict[str, Any] | None = None
) -> tuple[plt.Figure, npt.NDArray]:
    """
    Builds subplot figure and axes.

    Args:
        n_rows (int): Number of rows.
        n_cols (int): Number of columns.
        fig_params (Dict[str: Any]): Matplotlib subplot params.

    Returns
        plt.Figure: Matplotlib figure.
        plt.Axes: Matplotlib axes.
    """
    if fig_params is None:
        fig_params = {
            "constrained_layout": True,
            "dpi": 200
        }
    if "dpi" not in fig_params.keys():
        fig_params["dpi"] = 200
    if "constrained_layout" not in fig_params.keys():
        fig_params["constrained_layout"] = True

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        **fig_params
    )

    return fig, axes
