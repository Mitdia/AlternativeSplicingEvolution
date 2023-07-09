TYPE_SUBSTITUTION_PRICE = 1e+4
NEW_SPLICE_SITE_PRICE = 1e+3
NEW_SPLICING_VARIANT_PRICE = 1e+4
INTRON_TO_EXON_PRICE = 1e+6


DELETION_PRICE = 0
EXON_TO_INTRON = 1
START_EMERGENCE = 1
END_EMERGENCE = 1
START_TO_END = 1e+10
START_TO_EXON = 1
START_TO_INTRON = 1
END_TO_EXON = 1
END_TO_INTRON = 1
SHORTER_ISOFORM_PRICE = 1

SUBST_PRICES = {"x": {"x": 0, "-": DELETION_PRICE, "e": DELETION_PRICE, "S": DELETION_PRICE, "F": DELETION_PRICE, "=": SHORTER_ISOFORM_PRICE},
                "-": {"x": DELETION_PRICE, "-": 0, "e": EXON_TO_INTRON, "S": START_EMERGENCE, "F": END_EMERGENCE, "=": SHORTER_ISOFORM_PRICE},
                "e": {"x": DELETION_PRICE, "-": EXON_TO_INTRON, "e": 0, "S": START_EMERGENCE, "F": END_EMERGENCE, "=": SHORTER_ISOFORM_PRICE},
                "S": {"x": DELETION_PRICE, "-": START_TO_INTRON, "e": START_TO_EXON, "S": 0, "F": START_TO_END, "=": SHORTER_ISOFORM_PRICE},
                "F": {"x": DELETION_PRICE, "-": END_TO_INTRON, "e": END_TO_EXON, "S": START_TO_END, "F": 0, "=": SHORTER_ISOFORM_PRICE},
                "=": {"x": SHORTER_ISOFORM_PRICE, "-": SHORTER_ISOFORM_PRICE, "e": SHORTER_ISOFORM_PRICE, "S": SHORTER_ISOFORM_PRICE, "F": SHORTER_ISOFORM_PRICE, "=": 0}}


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
