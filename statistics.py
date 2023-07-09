import math
from matplotlib  import pyplot as plt


def percentage_of_coinciding_splice_sites(first_gene_structures, second_gene_structures):
    average_percentage_of_coinciding_donor_site = []
    average_percentage_of_coinciding_acceptor_site = []
    average_percentage_of_coinciding_splice_site = []
    for first_gene_donors, first_gene_acceptors in first_gene_structures:
        for second_gene_donors, second_gene_acceptors in second_gene_structures:
            number_of_donor_sites = len((set(first_gene_donors) | set(second_gene_donors)))
            number_of_coinciding_donor_sites = len((set(first_gene_donors) & set(second_gene_donors)))
            average_percentage_of_coinciding_donor_site.append(number_of_coinciding_donor_sites / number_of_donor_sites)
            number_of_acceptor_sites = len((set(first_gene_acceptors) | set(second_gene_acceptors)))
            number_of_coinciding_acceptor_sites = len((set(first_gene_acceptors) & set(second_gene_acceptors)))
            average_percentage_of_coinciding_acceptor_site.append(number_of_coinciding_acceptor_sites / number_of_acceptor_sites)
            number_of_splice_sites = number_of_acceptor_sites + number_of_donor_sites
            number_of_coinciding_splice_sites = number_of_coinciding_acceptor_sites + number_of_coinciding_donor_sites
            average_percentage_of_coinciding_splice_site.append(number_of_coinciding_splice_sites / number_of_splice_sites)
    return average_percentage_of_coinciding_donor_site, average_percentage_of_coinciding_acceptor_site, average_percentage_of_coinciding_splice_site


def distance_between_closest_splice_sites(first_gene_structures, second_gene_structures):
    first_gene_donor_sites = set()
    first_gene_acceptor_sites = set()
    second_gene_donor_sites = set()
    second_gene_acceptor_sites = set()

    for donors, acceptors in first_gene_structures:
        first_gene_acceptor_sites |= set(acceptors)
        first_gene_donor_sites |= set(donors)
    for donors, acceptors in second_gene_structures:
        second_gene_acceptor_sites |= set(acceptors)
        second_gene_donor_sites |= set(donors)

    donor_distances = []
    for first_gene_donor_site in first_gene_donor_sites:
        min_distance = math.inf
        for second_gene_donor_site in second_gene_donor_sites:
            distance = abs(second_gene_donor_site - first_gene_donor_site)
            if distance < min_distance:
                min_distance = distance
        donor_distances.append(min_distance)

    acceptor_distances = []
    for first_gene_acceptor_site in first_gene_acceptor_sites:
        min_distance = math.inf
        for second_gene_acceptor_site in second_gene_acceptor_sites:
            distance = abs(second_gene_acceptor_site - first_gene_acceptor_site)
            if distance < min_distance:
                min_distance = distance
        acceptor_distances.append(min_distance)

    return donor_distances, acceptor_distances


def visualize(closest_donor_distances, closest_acceptor_distances,
              coinciding_donors, coinciding_acceptors, coinciding_splice_sites):
    plt.ylabel("Number of splice sites")
    plt.xlabel("Distance between closest splice sites")
    plt.title(f"Histogram of distances between closest donor sites\n"
              f"(excluding purly coinciding "
              f"{closest_donor_distances.count(0)} of "
              f"{len(closest_donor_distances)})")
    plt.grid()
    plt.hist(closest_donor_distances, bins=300, range=(1, 300))
    plt.show()

    plt.ylabel("Number of splice sites")
    plt.xlabel("Distance between closest splice sites")
    plt.title(f"Histogram of distances between closest acceptor sites\n"
              f"(excluding purly coinciding "
              f"{closest_acceptor_distances.count(0)} of "
              f"{len(closest_acceptor_distances)})")
    plt.grid()
    plt.hist(closest_acceptor_distances, bins=300, range=(1, 300))
    plt.show()

    plt.ylabel("Number of splice sites")
    plt.xlabel("Distance between closest splice sites")
    plt.title(f"Histogram of distances between closest splice sites\n"
              f"(excluding purly coinciding "
              f"{(closest_donor_distances + closest_acceptor_distances).count(0)} of "
              f"{len(closest_donor_distances + closest_acceptor_distances)})")
    plt.grid()
    plt.hist(closest_donor_distances + closest_acceptor_distances, bins=300, range=(1, 300))
    plt.show()

    plt.xlabel("Percentage of coinciding donor sites")
    plt.ylabel("Number of pairs")
    plt.title("Number of pairs by percentage of coinciding donor sites")
    plt.grid()
    plt.hist(coinciding_donors, bins=100)
    plt.show()

    plt.xlabel("Percentage of coinciding acceptor sites")
    plt.ylabel("Number of pairs")
    plt.title("Number of pairs by percentage of coinciding acceptor sites")
    plt.grid()
    plt.hist(coinciding_acceptors, bins=100)
    plt.show()

    plt.xlabel("Percentage of coinciding splice sites")
    plt.ylabel("Number of pairs")
    plt.title("Number of pairs by percentage of coinciding splice sites")
    plt.grid()
    plt.hist(coinciding_splice_sites, bins=100)
    plt.show()
