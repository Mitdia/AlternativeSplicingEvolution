import re


def structure_to_coordinate_list(structure):
    exons = re.finditer(r"[eSFx]*e[eSFx]*", structure)
    acceptors = []
    donors = []
    for exon in exons:
        acceptors.append(exon.start())
        donors.append(exon.end())
    return donors, acceptors


def parse_file(filepath):
    with open(filepath, 'r') as file:

        file_contents = file.read()
        file_contents = file_contents.split()
        length = len(file_contents[1])
        at_structures = {}
        i = 2

        while file_contents[i].startswith(">AT"):
            structure = file_contents[i + 1]
            at_structures[structure] = (structure_to_coordinate_list(structure))
            i += 2

        i += 2
        bra_structures = {}
        while i < len(file_contents):
            structure = file_contents[i + 1]
            bra_structures[structure] = (structure_to_coordinate_list(structure))
            i += 2

        return at_structures, bra_structures
