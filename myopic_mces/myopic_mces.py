# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 17:16:05 2020
@author: seipp
"""
import time

import networkx as nx

from typing import Tuple

from myopic_mces.graph import construct_graph
from myopic_mces.MCES_ILP import MCES_ILP
from myopic_mces.filter_MCES import apply_filter


# noinspection PyPep8Naming
def MCES_core(ind: int, g1: nx.Graph, g2: nx.Graph, threshold: int, solver: str, solver_options=None, no_ilp_threshold=False, always_stronger_bound=True) -> Tuple[int, float, int]:
    if solver_options is None:
        solver_options = {}

    # filter out if distance is above the threshold
    if threshold == -1:
        res = MCES_ILP(g1, g2, threshold, solver, solver_options=solver_options, no_ilp_threshold=no_ilp_threshold)
        return ind, res[0], res[1]

    d, filter_id = apply_filter(g1, g2, threshold, always_stronger_bound=always_stronger_bound)

    if d > threshold:
        return ind, d, filter_id

    # calculate MCES
    res = MCES_ILP(g1, g2, threshold, solver, solver_options=solver_options, no_ilp_threshold=no_ilp_threshold)

    return ind, res[0], res[1]


# noinspection PyPep8Naming
# noinspection GrazieInspection
def MCES(ind: int, s1: str, s2: str, threshold: int, solver: str, solver_options=None, no_ilp_threshold=False, always_stronger_bound=True) -> Tuple[int, float, float, int]:
    """
    Calculates the distance between two molecules

    Parameters
    ----------
    ind : int
        index
    s1 : str
        SMILES of the first molecule
    s2 : str
        SMILES of the second molecule
    threshold : int
        Threshold for the comparison. Exact distance is only calculated if the distance is lower than the threshold.
        If set to -1 the exact distance is always calculated.
    solver : string
        ILP-solver used for solving MCES. Example:CPLEX_CMD
    solver_options : dict
        additional options to pass to solvers. Example: threads=1 for better multi-threaded performance
    no_ilp_threshold : bool
        if true, always return exact distance even if it is below the threshold (slower)
    always_stronger_bound : bool
        if true, always compute and use the second stronger bound

    Returns
    -------
    int
        index
    float
        Distance between the molecules
    float
        Time taken for the calculation
    int
        Type of Distance:
            1 : Exact Distance
            2 : Lower bound (if the exact distance is above the threshold; bound chosen dynamically)
            4 : Lower bound (second lower bound was used)
    """
    start = time.time()

    if solver_options is None:
        solver_options = {}

    # construct graph for both smiles.
    g1 = construct_graph(s1)
    g2 = construct_graph(s2)

    ind, res0, res1 = MCES_core(ind, g1, g2, threshold, solver, solver_options, no_ilp_threshold, always_stronger_bound)

    end = time.time()
    total_time = end - start

    return ind, res0, total_time, res1
