# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 17:16:05 2020
@author: seipp
"""
import time
import argparse
import multiprocessing

import networkx as nx

from typing import Tuple, List
from joblib import Parallel, delayed

from src.myopic_mces.graph import construct_graph
from src.myopic_mces.MCES_ILP import MCES_ILP
from src.myopic_mces.filter_MCES import apply_filter


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


def read_input_file(input_file_path: str) -> List[Tuple[int, str, str]]:
    with open(input_file_path, 'r') as f:
        inputs = []
        for line in f:
            args = line.split(',')
            inputs.append((int(args[0]), args[1], args[2]))
        return inputs


def write_output_file(output_file_path: str, results: List[Tuple[int, float, float, int]]):
    with open(output_file_path, 'w') as out:
        for i in results:
            out.write(f'{i[0]},{i[2]},{i[1]},{i[3]}\n')
    pass


def main_core(args):
    additional_mces_options = dict(no_ilp_threshold=args.no_ilp_threshold, solver_options=dict(), always_stronger_bound=not args.choose_bound_dynamically)

    if args.solver_onethreaded:
        additional_mces_options['solver_options']['threads'] = 1
    if args.solver_no_msg:
        additional_mces_options['solver_options']['msg'] = False

    inputs = read_input_file(args.input)

    num_jobs = multiprocessing.cpu_count() if args.num_jobs is None else args.num_jobs
    if num_jobs > 1:
        results = Parallel(n_jobs=num_jobs, verbose=5)(delayed(MCES)(ind, s1, s2, args.threshold, args.solver, **additional_mces_options) for ind, s1, s2 in inputs)
    else:
        results = [MCES(ind, s1, s2, args.threshold, args.solver, **additional_mces_options) for ind, s1, s2 in inputs]

    write_output_file(args.output, results)
    pass


# noinspection SpellCheckingInspection
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input", help="input file in the format: id,smiles1,smiles2")
    parser.add_argument("output", help="output file")
    parser.add_argument("--threshold", type=int, default=10,
                        action="store", help="threshold for the distance")
    parser.add_argument("--no_ilp_threshold", action="store_true",
                        help="(experimental) if set, do not add threshold as constraint to ILP, "
                        "resulting in longer runtimes and potential violations of the triangle equation")
    parser.add_argument("--choose_bound_dynamically", action="store_true",
                        help="if this is set, compute and use potentially weaker but faster lower bound if "
                        "already greater than the threshold. Otherwise (default), the strongest lower bound "
                        "is always computed and used")
    parser.add_argument("--solver", type=str, default="default",
                        action="store", help="Solver for the ILP. example:CPLEX_CMD")
    parser.add_argument("--solver_onethreaded", action="store_true",
                        help="limit ILP solver to one thread, resulting in faster "
                        "performance with parallel computations (not available for all solvers)")
    parser.add_argument("--solver_no_msg", action="store_true",
                        help="prevent solver from logging (not available for all solvers)")
    parser.add_argument("--num_jobs", type=int, help="Number of jobs; instances to run in parallel. "
                        "By default this is set to the number of (logical) CPU cores.")
    args = parser.parse_args()

    main_core(args)
    pass


if __name__ == '__main__':
    main()
