// Writing code for dynamic louvain community detection algorithm
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

double calculate_modularity_change(const vector<vector<pair<int, double>>>& adj, const vector<int>& community, int node, int old_community, int new_community, vector<double>& community_degree, vector<double>& community_internal) {
    double m = 0.0;
    for (const auto& neighbors : adj) {
        for (const auto& [_, w] : neighbors) {
            m += w;
        }
    }

    double delta_Q = 0.0;
    double k_i_in_old = 0.0, k_i_in_new = 0.0;
    double k_i = 0.0;
    double sum_tot_old = community_degree[old_community], sum_tot_new = community_degree[new_community];

    for (const auto& [v, w] : adj[node]) {
        k_i += w;
        if (community[v] == old_community) {
            k_i_in_old += w;
        }
        if (community[v] == new_community) {
            k_i_in_new += w;
        }
    }

    cout << 2.0 * (k_i_in_new - k_i_in_old) / m  << "This is internal changes\n";
    cout << "sum_tot_old: " << sum_tot_old << " sum_tot_new: " << sum_tot_new << " k_i: " << k_i << "\n";
    cout << 2.0 * ((k_i * (sum_tot_old - sum_tot_new - k_i)) / (m * m)) << "This is degree changes\n";
    delta_Q = 2.0 * (k_i_in_new - k_i_in_old) / m + 2.0 * ((k_i * (sum_tot_old - sum_tot_new - k_i)) / (m * m));
    return delta_Q;
}

void update_community_degree_internal(const vector<vector<pair<int, double>>>& adj, const vector<int>& community, vector<double>& community_degree, vector<double>& community_internal, int node_to_move, int old_community, int new_community) {
    int n = adj.size();

    for (const auto& [v, w] : adj[node_to_move]) {
        if(community[v] == old_community) {
            community_internal[old_community] -= w;
        }
        else if(community[v] == new_community) {
            community_internal[new_community] += w;
        }
    }

    for (const auto& [v, w] : adj[node_to_move]) {
        community_degree[old_community] -= w;
        community_degree[new_community] += w;
    }
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
        // Imagine using emplace back instead of push back
        adj[u].emplace_back(v, w);
        adj[v].emplace_back(u, w);
    }
    vector<vector<pair<int, double>>> original_adjacency_list = adj;
    int orig_n = n;
    vector<int> community(n);
    vector<double> community_degree(n, 0.0);
    vector<double> community_internal(n, 0.0);
    for (int i = 0; i < n; i++) {
        community[i] = i;
        for(auto & [v, w] : adj[i]){
            community_degree[i] += w;
        }
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
                double modularity_change = calculate_modularity_change(adj, community, node_to_move, community[node_to_move], community[target_community_node], community_degree, community_internal);
                cout << "Try Move " << node_to_move << " from community " << community[node_to_move] << " to " << community[target_community_node] << "\n";
                print_communities(community, partition_nodes);
                cout << "Current modularity: " << current_modularity << " " ;
                cout << "Possible new modularity: " << modularity_change + current_modularity << "\n";
                if(modularity_change > 0){
                    current_modularity += modularity_change;
                    update_community_degree_internal(adj, community, community_degree, community_internal, node_to_move, community[node_to_move], community[target_community_node]);
                    community[node_to_move] = community[target_community_node];
                    changed = true;
                }
            }
        }
    }

    changed_community = aggregation(adj, community, partition_nodes, n);
    cout << "Aggregating Graph\n";
    cout << "Communities:\n";
    print_communities(community, partition_nodes);

    }

    cout << "Final Communities and Adjacency List of Partitions after Static phase:\n";
    print_communities(community, partition_nodes);
    print_adjacency_list(adj);
    cout << "Final Modularity: " << calculate_modularity(adj, community) << endl;

    // Now do dynamic louvain
    // Note that the partition nodes data structure is not really used anywhere, it is just an information structure
    adj = original_adjacency_list;
    n = orig_n;
    community = vector<int>(n);
    for(int partition_index = 0; partition_index < partition_nodes.size(); partition_index++){
        vector<int>& node_list = partition_nodes[partition_index];
        for(int node : node_list){
            community[node] = partition_index;
        }
    }
        // Resetting partitions. I define partition as just a group of nodes that we keep together.
    partition_nodes.resize(n);
    for(int i = 0; i < n; i++) {
        partition_nodes[i] = {i};
    }

    int num_queries;
    // check if the cin fails, if so end the program
    if (!(cin >> num_queries)) {
        cerr << "Error reading number of queries." << endl;
        return 1;
    }
    for(int q = 0; q < num_queries; q++){
        // Adding u, v edge with weight w
            // Here w can be negative or postive
        int u, v;
        double w;
        cin >> u >> v >> w;
        adj[u].emplace_back(v, w);
        adj[v].emplace_back(u, w);

        if(community[u] == community[v] && w >= 0){
            // Do nothing
        } else if(community[u] != community[v] && w <= 0){
            // Do nothing
        } else {
            // Try to see if any movement is possible
            
            queue<int> active_nodes;
            active_nodes.push(u);
            active_nodes.push(v);
            double current_modularity = calculate_modularity(adj, community);
            while(!active_nodes.empty()){
                int node_to_move = active_nodes.front();
                active_nodes.pop();
                int initial_community = community[node_to_move];
                for(int target_community_node = 0; target_community_node < n; target_community_node++){
                    if(community[node_to_move] == community[target_community_node])
                        continue;
                    community[node_to_move] = community[target_community_node];
                    double new_modularity = calculate_modularity(adj, community);
                    cout << "Possible new Community and modularity after dynamic change:\n";
                    print_communities(community, partition_nodes);
                    cout << new_modularity << "\n";
                    if(new_modularity > current_modularity){
                        current_modularity = new_modularity;
                        initial_community = community[node_to_move];
                        // Add all neighbors to active nodes
                        for(auto& [neighbor, _] : adj[node_to_move]){
                            active_nodes.push(neighbor);
                        }
                    } else {    
                        community[node_to_move] = initial_community;
                    }
                }
            }
        }

        cout << "Communities after handling query:\n";
        print_communities(community, partition_nodes);
    }


    changed_community = aggregation(adj, community, partition_nodes, n);
    cout << "Aggregating Graph\n";
    cout << "Communities:\n";
    print_communities(community, partition_nodes);
    
    cout << "Final Communities and Adjacency List of Partitions:\n";
    print_communities(community, partition_nodes);
    print_adjacency_list(adj);
    cout << "Final Modularity: " << calculate_modularity(adj, community) << endl;
    
}