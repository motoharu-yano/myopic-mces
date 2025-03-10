# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 17:17:41 2020

@author: seipp
"""
import pulp
import networkx as nx

from typing import Tuple


# noinspection PyPep8Naming
# noinspection GrazieInspection
def MCES_ILP(g1: nx.Graph, g2: nx.Graph, threshold, solver, solver_options=None, no_ilp_threshold=False) -> Tuple[float, int]:
    """
     Calculates the exact distance between two molecules using an ILP

     Parameters
     ----------
     g1 : networkx.classes.graph.Graph
         Graph representing the first molecule.
     g2 : networkx.classes.graph.Graph
         Graph representing the second molecule.
     threshold : int
         Threshold for the comparison. Exact distance is only calculated if the distance is lower than the threshold.
     solver: string
         ILP-solver used for solving MCES. Example:GUROBI_CMD
     solver_options: dict
         additional options to pass to solvers. Example: threads=1, msg=False for better multi-threaded performance
     no_ilp_threshold: bool
         if true, always return exact distance even if it is below the threshold (slower)

     Returns
     -------
     float
         Distance between the molecules
     int
         Type of Distance:
             1 : Exact Distance
             2 : Lower bound (If the exact distance is above the threshold)

    """
    if solver_options is None:
        solver_options = {}

    ILP = pulp.LpProblem("MCES", pulp.LpMinimize)

    # Variables for nodepairs
    nodepairs = []
    for i in g1.nodes:
        for j in g2.nodes:
            if g1.nodes[i]["atom"] == g2.nodes[j]["atom"]:
                nodepairs.append(tuple([i, j]))
    y = pulp.LpVariable.dicts('nodepairs', nodepairs,
                              lowBound=0,
                              upBound=1,
                              cat=pulp.LpInteger)

    # variables for edgepairs and weight
    edgepairs = []
    w = {}

    for i in g1.edges:
        for j in g2.edges:
            if (g1.nodes[i[0]]["atom"] == g2.nodes[j[0]]["atom"] and g1.nodes[i[1]]["atom"] == g2.nodes[j[1]]["atom"]) or (g1.nodes[i[1]]["atom"] == g2.nodes[j[0]]["atom"] and g1.nodes[i[0]]["atom"] == g2.nodes[j[1]]["atom"]):
                edgepairs.append(tuple([i, j]))
                w[tuple([i, j])] = max(g1[i[0]][i[1]]["weight"], g2[j[0]][j[1]]["weight"]) - min(g1[i[0]][i[1]]["weight"], g2[j[0]][j[1]]["weight"])

    # variables for not mapping an edge
    for i in g1.edges:
        edgepairs.append(tuple([i, -1]))
        w[tuple([i, -1])] = g1[i[0]][i[1]]["weight"]

    for j in g2.edges:
        edgepairs.append(tuple([-1, j]))
        w[tuple([-1, j])] = g2[j[0]][j[1]]["weight"]

    c = pulp.LpVariable.dicts('edgepairs', edgepairs,
                              lowBound=0,
                              upBound=1,
                              cat=pulp.LpInteger)

    # objective function
    ILP += pulp.lpSum([w[i]*c[i] for i in edgepairs])

    # Every node in G1 can only be mapped to at most one in G2
    for i in g1.nodes:
        h = []
        for j in g2.nodes:
            if g1.nodes[i]["atom"] == g2.nodes[j]["atom"]:
                h.append(tuple([i, j]))
        ILP += pulp.lpSum([y[k] for k in h]) <= 1

    # Every node in G1 can only be mapped to at most one in G1
    for i in g2.nodes:
        h = []
        for j in g1.nodes:
            if g1.nodes[j]["atom"] == g2.nodes[i]["atom"]:
                h.append(tuple([j, i]))
        ILP += pulp.lpSum([y[k] for k in h]) <= 1

    # Every edge in G1 has to be mapped to an edge in G2 or the variable for not mapping has to be 1
    for i in g1.edges:
        ls = []
        # rs = []
        for j in g2.edges:
            if (g1.nodes[i[0]]["atom"] == g2.nodes[j[0]]["atom"] and g1.nodes[i[1]]["atom"] == g2.nodes[j[1]]["atom"]) or (g1.nodes[i[1]]["atom"] == g2.nodes[j[0]]["atom"] and g1.nodes[i[0]]["atom"] == g2.nodes[j[1]]["atom"]):
                ls.append(tuple([i, j]))
        ILP += pulp.lpSum([c[k] for k in ls])+c[tuple([i, -1])] == 1

    # Every edge in G2 has to be mapped to an edge in G1 or the variable for not mapping has to be 1
    for i in g2.edges:
        ls = []
        # rs = []
        for j in g1.edges:
            if (g1.nodes[j[0]]["atom"] == g2.nodes[i[0]]["atom"] and g1.nodes[j[1]]["atom"] == g2.nodes[i[1]]["atom"]) or (g1.nodes[j[1]]["atom"] == g2.nodes[i[0]]["atom"] and g1.nodes[j[0]]["atom"] == g2.nodes[i[1]]["atom"]):
                ls.append(tuple([j, i]))
        ILP += pulp.lpSum([c[k] for k in ls])+c[tuple([-1, i])] == 1

    # The mapping of the edges has to match the mapping of the nodes
    for i in g1.nodes:
        for j in g2.edges:
            ls = []
            for k in g1.neighbors(i):
                if tuple([tuple([i, k]), j]) in c:
                    ls.append(tuple([tuple([i, k]), j]))
                else:
                    if tuple([tuple([k, i]), j]) in c:
                        ls.append(tuple([tuple([k, i]), j]))
            rs = []
            if g1.nodes[i]["atom"] == g2.nodes[j[0]]["atom"]:
                rs.append(tuple([i, j[0]]))
            if g1.nodes[i]["atom"] == g2.nodes[j[1]]["atom"]:
                rs.append(tuple([i, j[1]]))
            ILP += pulp.lpSum([c[k] for k in ls]) <= pulp.lpSum([y[k] for k in rs])

    for i in g2.nodes:
        for j in g1.edges:
            ls = []
            for k in g2.neighbors(i):
                if tuple([j, tuple([i, k])]) in c:
                    ls.append(tuple([j, tuple([i, k])]))
                else:
                    if tuple([j, tuple([k, i])]) in c:
                        ls.append(tuple([j, tuple([k, i])]))
            rs = []
            if g2.nodes[i]["atom"] == g1.nodes[j[0]]["atom"]:
                rs.append(tuple([j[0], i]))
            if g2.nodes[i]["atom"] == g1.nodes[j[1]]["atom"]:
                rs.append(tuple([j[1], i]))
            ILP += pulp.lpSum([c[k] for k in ls]) <= pulp.lpSum(y[k] for k in rs)

    # constraint for the threshold
    if threshold != -1 and not no_ilp_threshold:
        ILP += pulp.lpSum([w[i]*c[i] for i in edgepairs]) <= threshold

    # solve the ILP
    if solver == "default":
        ILP.solve()
    else:
        sol = pulp.getSolver(solver, **solver_options)
        ILP.solve(sol)
    if ILP.status == 1:
        return float(ILP.objective.value()), 1
    else:
        return threshold, 2
