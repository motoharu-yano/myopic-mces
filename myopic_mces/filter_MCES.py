# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 17:59:05 2020

@author: seipp
"""

import networkx as nx

from typing import Tuple


def filter1(g1: nx.Graph, g2: nx.Graph) -> float:
    """
     Finds a lower bound for the distance based on degree

     Parameters
     ----------
     g1 : networkx.classes.graph.Graph
         Graph representing the first molecule.
     g2 : networkx.classes.graph.Graph
         Graph representing the second molecule.

     Returns
     -------
     float
         Lower bound for the distance between the molecules

    """
    # Find all occurring atom types and partition by type
    atom_types1 = []
    for i in g1.nodes:
        if g1.nodes[i]["atom"] not in atom_types1:
            atom_types1.append(g1.nodes[i]["atom"])

    type_map1 = {}
    for i in atom_types1:
        type_map1[i] = list(filter(lambda x: i == g1.nodes[x]["atom"], g1.nodes))

    atom_types2 = []
    for i in g2.nodes:
        if g2.nodes[i]["atom"] not in atom_types2:
            atom_types2.append(g2.nodes[i]["atom"])

    type_map2 = {}
    for i in atom_types2:
        type_map2[i] = list(filter(lambda x: i == g2.nodes[x]["atom"], g2.nodes))

    # calculate lower bound
    difference = 0

    # Every atom type is done seperately
    for i in atom_types1:
        if i in atom_types2:
            # number of nodes that can be mapped
            n = min(len(type_map1[i]), len(type_map2[i]))

            # sort by degree
            degreelist1 = sorted(type_map1[i], key=lambda x: sum([g1[x][j]["weight"] for j in g1.neighbors(x)]), reverse=True)
            degreelist2 = sorted(type_map2[i], key=lambda x: sum([g2[x][j]["weight"] for j in g2.neighbors(x)]), reverse=True)

            # map in order of sorted lists
            for j in range(n):
                deg1 = sum([g1[degreelist1[j]][k]["weight"] for k in g1.neighbors(degreelist1[j])])
                deg2 = sum([g2[degreelist2[j]][k]["weight"] for k in g2.neighbors(degreelist2[j])])
                difference += abs(deg1-deg2)
            # nodes that are not mapped
            if len(degreelist1) > n:
                for j in range(n, len(degreelist1)):
                    difference += sum([g1[degreelist1[j]][k]["weight"] for k in g1.neighbors(degreelist1[j])])
            if len(degreelist2) > n:
                for j in range(n, len(degreelist2)):
                    difference += sum([g2[degreelist2[j]][k]["weight"] for k in g2.neighbors(degreelist2[j])])
        # atom type only in one of the graphs
        else:
            for j in type_map1[i]:
                difference += sum([g1[j][k]["weight"] for k in g1.neighbors(j)])
    for i in atom_types2:
        if i not in atom_types1:
            for j in type_map2[i]:
                difference += sum([g2[j][k]["weight"] for k in g2.neighbors(j)])
    return difference/2


def get_cost(g1: nx.Graph, g2: nx.Graph, i: int, j: int) -> float:
    """
     Calculates the cost for mapping node i to j based on neighborhood

     Parameters
     ----------
     g1 : networkx.classes.graph.Graph
         Graph representing the first molecule.
     g2 : networkx.classes.graph.Graph
         Graph representing the second molecule.
     i : int
         Node of G1
     j : int
         Node of G2

     Returns
     -------
     float
         Cost of mapping i to j

    """
    # Find all occuring atom types in neighborhood
    atom_types1 = []
    for k in g1.neighbors(i):
        if g1.nodes[k]["atom"] not in atom_types1:
            atom_types1.append(g1.nodes[k]["atom"])
    type_map1 = {}
    for k in atom_types1:
        type_map1[k] = list(filter(lambda x: k == g1.nodes[x]["atom"], g1.neighbors(i)))

    atom_types2 = []
    for k in g2.neighbors(j):
        if g2.nodes[k]["atom"] not in atom_types2:
            atom_types2.append(g2.nodes[k]["atom"])

    type_map2 = {}
    for k in atom_types2:
        type_map2[k] = list(filter(lambda x: k == g2.nodes[x]["atom"], g2.neighbors(j)))

    # calculate cost
    difference = 0.0

    # Every atom type is handled seperately
    for k in atom_types1:
        if k in atom_types2:
            n = min(len(type_map1[k]), len(type_map2[k]))
            # sort by incident edges by weight
            edgelist1 = sorted(type_map1[k], key=lambda x: g1[i][x]["weight"], reverse=True)
            edgelist2 = sorted(type_map2[k], key=lambda x: g2[j][x]["weight"], reverse=True)

            # map in order of sorted lists
            for l in range(n):
                difference += (max(g1[i][edgelist1[l]]["weight"], g2[j][edgelist2[l]]["weight"]) - min(g1[i][edgelist1[l]]["weight"], g2[j][edgelist2[l]]["weight"])) / 2

            # cost for not mapped edges
            if len(edgelist1) > n:
                for l in range(n, len(edgelist1)):
                    difference += g1[i][edgelist1[l]]["weight"] / 2

            if len(edgelist2) > n:
                for l in range(n, len(edgelist2)):
                    difference += g2[j][edgelist2[l]]["weight"] / 2
        else:
            for l in type_map1[k]:
                difference += g1[i][l]["weight"] / 2

    for k in atom_types2:
        if k not in atom_types1:
            for l in type_map2[k]:
                difference += g2[j][l]["weight"] / 2
    return difference


def filter2(g1: nx.Graph, g2: nx.Graph) -> float:
    """
     Finds a lower bound for the distance based on neighborhood

     Parameters
     ----------
     g1 : networkx.classes.graph.Graph
         Graph representing the first molecule.
     g2 : networkx.classes.graph.Graph
         Graph representing the second molecule.

     Returns
     -------
     float
         Lower bound for the distance between the molecules

    """
    # Find all occuring atom types
    atom_types1 = []
    for i in g1.nodes:
        if g1.nodes[i]["atom"] not in atom_types1:
            atom_types1.append(g1.nodes[i]["atom"])

    atom_types2 = []
    for i in g2.nodes:
        if g2.nodes[i]["atom"] not in atom_types2:
            atom_types2.append(g2.nodes[i]["atom"])

    atom_types = atom_types1

    for i in atom_types2:
        if i not in atom_types:
            atom_types.append(i)

    # calculate distance
    res = 0.0

    # handle every atom type seperately
    for i in atom_types:
        # filter by atom type
        nodes1 = list(filter(lambda x: i == g1.nodes[x]["atom"], g1.nodes))
        nodes2 = list(filter(lambda x: i == g2.nodes[x]["atom"], g2.nodes))

        # Create new graph for and solve minimum weight full matching
        G = nx.Graph()

        # Add node for every node of type i in G1 and G2
        for j in nodes1:
            G.add_node(tuple([1, j]))

        for j in nodes2:
            G.add_node(tuple([2, j]))

        # Add edges between all nodes of G1 and G2
        for j in nodes1:
            for k in nodes2:
                if g1.nodes[j]["atom"] == g2.nodes[k]["atom"]:
                    G.add_edge(tuple([1, j]), tuple([2, k]), weight=get_cost(g1, g2, j, k))

        # Add nodes if one graph has more nodes of type i than the other
        if len(nodes1) < len(nodes2):
            diff = len(nodes2)-len(nodes1)
            for j in range(1, diff+1):
                G.add_node(tuple([1, -j]))
                for k in nodes2:
                    G.add_edge(tuple([1, -j]), tuple([2, k]), weight=sum([g2[l][k]["weight"] for l in g2.neighbors(k)]) / 2)

        if len(nodes2) < len(nodes1):
            diff = len(nodes1)-len(nodes2)
            for j in range(1, diff+1):
                G.add_node(tuple([2, -j]))
                for k in nodes1:
                    G.add_edge(tuple([1, k]), tuple([2, -j]), weight=sum([g1[l][k]["weight"] for l in g1.neighbors(k)]) / 2)

        # Solve minimum weight full matching
        h = nx.bipartite.minimum_weight_full_matching(G)

        # Add weight of the matching
        for k in h:
            if k[0] == 1:
                res = res+G[k][h[k]]["weight"]

    return res


def apply_filter(g1: nx.Graph, g2: nx.Graph, threshold: int, always_stronger_bound=True) -> Tuple[float, int]:
    """
     Finds a lower bound for the distance

     Parameters
     ----------
     g1 : networkx.classes.graph.Graph
         Graph representing the first molecule.
     g2 : networkx.classes.graph.Graph
         Graph representing the second molecule.
     threshold : int
         Threshold for the comparison. We want to find a lower bound that is higher than the threshold
     always_stronger_bound : bool
         if true, always compute and use the second stronger bound

     Returns
     -------
     float
         Lower bound for the distance between the molecules
     int
         Which lower bound was chosen: 2 - depending on threshold, 4 - second lower bound

    """
    if always_stronger_bound:
        d = filter2(g1, g2)
        return d, 4
    else:
        # calculate first lower bound
        d = filter1(g1, g2)
        # if below threshold calculate second lower bound
        if d <= threshold:
            d = filter2(g1, g2)
            if d <= threshold:
                return d, 2

        return d, 2
