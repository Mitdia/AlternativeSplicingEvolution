TYPE_SUBSTITUTION_PRICE = 1e+4
NEW_SPLICE_SITE_PRICE = 1e+3
NEW_SPLICING_VARIANT_PRICE = 1e+4
INTRON_TO_EXON_PRICE = 1e+6


def node_substitution_price(first_node, second_node):
    if first_node["type"] == second_node["type"]:
        return (first_node["coordinate"] - second_node["coordinate"]) ** 2
    else:
        return TYPE_SUBSTITUTION_PRICE


def node_mathing(first_node, second_node):
    if first_node["type"] == second_node["type"] and first_node["coordinate"] == second_node["coordinate"]:
        return True
    return False


def edge_substitution_price(first_edge, second_edge):
    if first_edge["type"] == second_edge["type"]:
        return 0
    return INTRON_TO_EXON_PRICE


def edge_match(first_edge, second_edge):
    if first_edge["type"] == second_edge["type"]:
        return True
    return False
