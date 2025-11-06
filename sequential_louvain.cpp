// Writing code for sequential louvain community detection algorithm
#include <bits/stdc++.h>
using namespace std;

double calculate_modularity(const vector<vector<pair<int, double>>>& adj, const vector<int>& community) {
    double m = 0.0;
    for (const auto& neighbors : adj) {
        for (const auto& [_, w] : neighbors) {
            m += w;
        }
    }

    double Q = 0.0;
    map<int, double> community_degree;
    map<int, double> community_internal;

    for (int u = 0; u < adj.size(); u++) {
        int comm_u = community[u];
        double degree_u = 0.0;
        for (const auto& [v, w] : adj[u]) {
            degree_u += w;
            if (community[v] == comm_u) {
                community_internal[comm_u] += w;
            }
        }
        community_degree[comm_u] += degree_u;
    }

    for (const auto& [comm, degree] : community_degree) {
        double internal = community_internal[comm];
        Q += (internal / (m)) - ((degree / (m)) * (degree / (m)));
    }

    return Q;
}

void print_communities(const vector<int>& community, const vector<vector<int>>& partition_nodes) {
    for (int i = 0; i < community.size(); i++) {
        cout << "Partition " << i << " in community " << community[i] << ": Nodes in partition - ";
        for(int node : partition_nodes[i]){
            cout << node << " ";
        }
        cout << "\n";
    }
}

void print_adjacency_list(const vector<vector<pair<int, double>>>& adj) {
    for (int u = 0; u < adj.size(); u++) {
        cout << "Node " << u << ": ";
        for (const auto& [v, w] : adj[u]) {
            cout << "(" << v << ", " << w << ") ";
        }
        cout << "\n";
    }
}

bool aggregation(vector<vector<pair<int, double>>>& adj, vector<int>& community, vector<vector<int>>& partition_nodes, int& n){
    map<int, int> community_map;
    int new_index = 0;
    bool community_movement = false; // If each community comes up only once, then there was no movement in the partitions/nodes
    for(int i = 0; i < n; i++){
        if(community_map.find(community[i]) == community_map.end()){
            community_map.insert({community[i], new_index});
            new_index++;
        } else {
            community_movement = true;
        }
        community[i] = community_map[community[i]];
    }
    // Created the communities from 0 to new_index - 1

    vector<map<int, double>> new_adj_map(new_index);
    for(int u = 0; u < n; u++){
        int new_u = community[u];
        for(const auto& [v, w] : adj[u]){
            int new_v = community[v];
            new_adj_map[new_u][new_v] += w;
        }
    }

    vector<vector<pair<int, double>>> new_adj(new_index);
    for(int u = 0; u < new_index; u++){
        for(const auto& [v, w] : new_adj_map[u]){
            new_adj[u].push_back({v, w});
        }
    }
    // Created new adjacency list
    
    vector<vector<int>> new_partition_nodes(new_index);
    for(int i = 0; i < n; i++){
        new_partition_nodes[community[i]].insert(new_partition_nodes[community[i]].end(), partition_nodes[i].begin(), partition_nodes[i].end());
    }
    // Created new partition nodes

    vector<int> new_community(new_index);
    for(int i = 0; i < new_index; i++){
        new_community[i] = i;
    }

    adj = new_adj;
    partition_nodes = new_partition_nodes;
    community = new_community;
    n = community.size();
    return community_movement;
}


int main(){
    int n, m;
    cin >> n >> m;
    vector<vector<pair<int, double>>> adj(n);
    vector<vector<int>> partition_nodes(n);
    for (int i = 0; i < m; i++) {
        int u, v;
        double w;
        cin >> u >> v >> w;
        // w = 1.0; // Assuming unweighted graph for simplicity
        // Imagine using emplace back instead of push back
        adj[u].emplace_back(v, w);
        adj[v].emplace_back(u, w);
    }
    vector<int> community(n);
    for (int i = 0; i < n; i++) {
        community[i] = i;
        partition_nodes[i].push_back(i);
    }

    bool changed_community = true;
    while(changed_community){
    changed_community = false; 

    bool changed = true;
    double current_modularity = calculate_modularity(adj, community);
    while(changed){
        // int node_moved = 0;
        // int community_target = 0;
        // Can just do dynamically

        changed = false;
        for(int node_to_move = 0; node_to_move < n; node_to_move++){
            int initial_community = community[node_to_move];
            for(int target_community_node = 0; target_community_node < n; target_community_node++){
                if(community[node_to_move] == community[target_community_node])
                    continue;
                community[node_to_move] = community[target_community_node];
                double new_modularity = calculate_modularity(adj, community);
                cout << "Possible new Community and modularity:\n";
                cout << "Move node " << node_to_move << " to community " << community[target_community_node] << "\n";
                // print_communities(community, partition_nodes);
                cout << new_modularity << "\n";
                if(new_modularity > current_modularity){
                    current_modularity = new_modularity;
                    initial_community = community[node_to_move];
                    changed = true;
                } else {    
                    community[node_to_move] = initial_community;
                }
            }
        }
    }

    changed_community = aggregation(adj, community, partition_nodes, n);
    cout << "Aggregating Graph\n";
    cout << "Communities:\n";
    print_communities(community, partition_nodes);

    }

    cout << "Final Communities and Adjacency List of Partitions:\n";
    print_communities(community, partition_nodes);
    print_adjacency_list(adj);

    cout << "Final Modularity: " << calculate_modularity(adj, community) << endl;
}