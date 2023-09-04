import numpy as np


def edge_length(node_i, node_j):
    return abs(node_i[0] - node_j[0]) + abs(node_i[1] - node_j[1])


class Node:
    def __init__(self, i, data):
        self.data = data
        self.i = i
        self.neighbors = []
        self.visited = False


def dfs(node):
    node.visited = True
    curr_comp = [node.data]
    for neighbor in node.neighbors:
        if not neighbor.visited:
            curr_comp.extend(dfs(neighbor))
    return curr_comp


def find_mean(component):
    # print(component)
    mean_coords = np.mean(component, axis=0)
    distances = np.sum(np.abs(np.array(component) - mean_coords[np.newaxis]), axis=1)
    min_idx = np.argmin(distances)
    return component[min_idx]


def reduce_clusters(all_nodes):
    node_list = [Node(i, x) for i, x in enumerate(all_nodes)]
    for i in range(len(all_nodes)):
        for j in range(i + 1, len(all_nodes)):
            if edge_length(all_nodes[i], all_nodes[j]) < 2:
                node_list[i].neighbors.append(node_list[j])
                node_list[j].neighbors.append(node_list[i])
    all_means = []
    for i in range(len(all_nodes)):
        if not node_list[i].visited:
            curr_component = dfs(node_list[i])
            curr_mean = find_mean(curr_component)
            all_means.append(curr_mean)
    return all_means
