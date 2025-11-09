import networkx as nx
import random

def create_weighted_random_graph_and_output_edgelist(num_nodes, probability_of_edge, min_weight, max_weight):
    # Create an Erdos-Renyi random graph
    # gnp_random_graph also can be used
    G = nx.erdos_renyi_graph(num_nodes, probability_of_edge)

    # Assign random weights to the edges
    for u, v in G.edges():
        G[u][v]['weight'] = random.randint(min_weight, max_weight)

    # Prepare the output string
    output_lines = []
    output_lines.append(f"{G.number_of_nodes()} {G.number_of_edges()}")

    for u, v, data in G.edges(data=True):
        output_lines.append(f"{u} {v} {data['weight']}")

    return "\n".join(output_lines)

# Example usage:
num_nodes = 100
probability_of_edge = 0.04
min_weight = 1
max_weight = 10

for i in range(5):
    edge_list_output = create_weighted_random_graph_and_output_edgelist(
        num_nodes, probability_of_edge, min_weight, max_weight
    )
    with open(f"./nxgraphs/nxgraph{i}", "w") as f:
        f.write(edge_list_output)

num_nodes = 1000
probability_of_edge = 0.02
min_weight = 1
max_weight = 10

for i in range(5, 10):
    edge_list_output = create_weighted_random_graph_and_output_edgelist(
        num_nodes, probability_of_edge, min_weight, max_weight
    )
    with open(f"./nxgraphs/nxgraph{i}", "w") as f:
        f.write(edge_list_output)