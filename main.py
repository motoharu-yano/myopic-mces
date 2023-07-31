# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 17:16:05 2020
@author: seipp
"""
import argparse
import multiprocessing

from typing import Tuple, List
from joblib import Parallel, delayed

from myopic_mces import myopic_mces as mces


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
        results = Parallel(n_jobs=num_jobs, verbose=5)(delayed(mces.MCES)(ind, s1, s2, args.threshold, args.solver, **additional_mces_options) for ind, s1, s2 in inputs)
    else:
        results = [mces.MCES(ind, s1, s2, args.threshold, args.solver, **additional_mces_options) for ind, s1, s2 in inputs]

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
