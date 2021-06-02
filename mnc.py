import scipy.sparse as sps
import numpy as np

def get_counterpart(alignment_matrix):
    counterpart_dict = {}

    if not sps.issparse(alignment_matrix):
        sorted_indices = np.argsort(alignment_matrix)

    n_nodes = alignment_matrix.shape[0]
    for node_index in range(n_nodes):

        if sps.issparse(alignment_matrix):
            row, possible_alignments, possible_values = sps.find(alignment_matrix[node_index])
            node_sorted_indices = possible_alignments[possible_values.argsort()]
        else:
            node_sorted_indices = sorted_indices[node_index]
        counterpart = node_sorted_indices[-1]
        counterpart_dict[node_index] = counterpart
    return counterpart_dict


def score_MNC(alignment_matrix, adj1, adj2):
    mnc = 0
    if sps.issparse(alignment_matrix): alignment_matrix = alignment_matrix.toarray()
    if sps.issparse(adj1): adj1 = adj1.toarray()
    if sps.issparse(adj2): adj2 = adj2.toarray()
    counter_dict = get_counterpart(alignment_matrix)
    node_num = alignment_matrix.shape[0]

    for i in range(node_num):
        a = np.array(adj1[i, :])
        one_hop_neighbor = np.flatnonzero(a)
        b = np.array(adj2[counter_dict[i], :])
        # neighbor of counterpart
        new_one_hop_neighbor = np.flatnonzero(b)

        one_hop_neighbor_counter = []

        for count in one_hop_neighbor:
            one_hop_neighbor_counter.append(counter_dict[count])

        num_stable_neighbor = np.intersect1d(new_one_hop_neighbor, np.array(one_hop_neighbor_counter)).shape[0]
        union_align = np.union1d(new_one_hop_neighbor, np.array(one_hop_neighbor_counter)).shape[0]

        sim = float(num_stable_neighbor) / union_align
        mnc += sim

    mnc /= node_num
    return mnc
