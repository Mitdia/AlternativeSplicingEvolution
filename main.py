import os
import networkx as nx
from matplotlib import pyplot as plt
from data_parser import structure_to_coordinate_list, graph_from_splice_sites
from price_functions import node_substitution_price, edge_substitution_price
from price_functions import NEW_SPLICING_VARIANT_PRICE, NEW_SPLICE_SITE_PRICE


directory = "Alignment_structure"


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
