import os
import networkx as nx
from matplotlib import pyplot as plt
from data_parser import parse_file
from price_functions import node_substitution_price, edge_substitution_price
from price_functions import NEW_SPLICING_VARIANT_PRICE, NEW_SPLICE_SITE_PRICE


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


directory = "Alignment_structure"


for filename in os.listdir(directory):
    filepath = os.path.join(directory, filename)
    if os.path.isfile(filepath):
        at_structures, bra_structures = parse_file(filepath)
        length = len(list(at_structures.keys())[0])
        at_graph = graph_from_splice_sites(at_structures.values(), length)
        bra_graph = graph_from_splice_sites(bra_structures.values(), length)
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
