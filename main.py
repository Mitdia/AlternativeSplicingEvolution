import os
import re
import networkx as nx
from matplotlib import pyplot as plt


TYPE_SUBSTITUTION_PRICE = 1e+4
NEW_SPLICE_SITE_PRICE = 1e+3
NEW_SPLICING_VARIANT_PRICE = 1e+4
INTRON_TO_EXON_PRICE = 1e+6


directory = "Alignment_structure"
data = []


def structure_to_coordinate_list(structure):
    exons = re.finditer(r"[eSFx]*e[eSFx]*", structure)
    acceptors = []
    donors = []
    for exon in exons:
        acceptors.append(exon.start())
        donors.append(exon.end())
    return donors, acceptors


def graph_from_splice_sites(structures, length):
    graph = nx.DiGraph()
    seen_sites = {}
    node_index = 1
    graph.add_node(0, coordinate=-1, type="root")
    graph.add_node(-1, coordinate=length + 1, type="end")
    for donors, acceptors in structures:

        if acceptors[0] not in seen_sites:
            graph.add_node(node_index, coordinate=acceptors[0], type="acceptor")
            seen_sites[acceptors[0]] = node_index
            node_index += 1
        graph.add_edge(0, seen_sites[acceptors[0]], type="intron")

        if donors[-1] not in seen_sites:
            graph.add_node(node_index, coordinate=donors[-1], type="donor")
            seen_sites[donors[-1]] = node_index
            node_index += 1
        graph.add_edge(seen_sites[donors[-1]], -1, type="intron")

        for i in range(1, len(donors)):

            if donors[i - 1] not in seen_sites:
                graph.add_node(node_index, coordinate=donors[i - 1], type="donor")
                seen_sites[donors[i - 1]] = node_index
                node_index += 1

            if acceptors[i] not in seen_sites:
                graph.add_node(node_index, coordinate=acceptors[i], type="acceptor")
                seen_sites[acceptors[i]] = node_index
                node_index += 1

            graph.add_edge(seen_sites[acceptors[i - 1]], seen_sites[donors[i - 1]], type="exon")
            graph.add_edge(seen_sites[donors[i - 1]], seen_sites[acceptors[i]], type="intron")

        graph.add_edge(seen_sites[acceptors[-1]], seen_sites[donors[-1]], type="exon")

    return graph


def node_substitution_price(first_node, second_node):
    if first_node["type"] == second_node["type"]:
        return (first_node["coordinate"] - second_node["coordinate"]) ** 2
    else:
        return TYPE_SUBSTITUTION_PRICE


def edge_substitution_price(first_edge, second_edge):
    if first_edge["type"] == second_edge["type"]:
        return 0
    return INTRON_TO_EXON_PRICE


def edge_match(first_edge, second_edge):
    if first_edge["type"] == second_edge["type"]:
        return True
    return False


def node_math(first_node, second_node):
    if first_node["type"] == second_node["type"] and first_node["coordinate"] == second_node["coordinate"]:
        return True
    return False


for filename in os.listdir(directory):
    filepath = os.path.join(directory, filename)
    if os.path.isfile(filepath):
        with open(filepath, 'r') as file:
            file_contents = file.read()
            file_contents = file_contents.split()
            length = len(file_contents[1])
            at_structures = []
            i = 2
            while file_contents[i].startswith(">AT"):
                at_structures.append(structure_to_coordinate_list(file_contents[i + 1]))
                i += 2
            at_graph = graph_from_splice_sites(at_structures, length)

            i += 2
            bra_structures = []
            while i < len(file_contents):
                bra_structures.append(structure_to_coordinate_list(file_contents[i + 1]))
                i += 2
            bra_graph = graph_from_splice_sites(bra_structures, length)
            nx.draw_networkx(at_graph)
            plt.show()
            nx.draw_networkx(bra_graph)
            plt.show()
            print(nx.graph_edit_distance(at_graph, bra_graph,
                                         node_subst_cost=node_substitution_price,
                                         node_ins_cost=lambda _: NEW_SPLICE_SITE_PRICE,
                                         node_del_cost=lambda _: NEW_SPLICE_SITE_PRICE,
                                         edge_subst_cost=edge_substitution_price,
                                         edge_ins_cost=lambda _: NEW_SPLICING_VARIANT_PRICE,
                                         edge_del_cost=lambda _: NEW_SPLICING_VARIANT_PRICE,
                                         roots=(0, 0),
                                         timeout=10))
