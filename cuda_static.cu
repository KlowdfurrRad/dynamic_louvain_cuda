// Code for cuda static louvain
#include <vector>
using namespace std;

#include <iostream>
#include <vector>
#include <cuda_runtime.h>
#include <cassert>

#include <thrust/device_vector.h> 
#include <thrust/transform.h> 
#include <thrust/sequence.h> 
#include <thrust/copy.h> 
#include <thrust/fill.h> 
#include <thrust/replace.h> 
#include <thrust/functional.h> 

#include <cooperative_groups.h>
namespace cg = cooperative_groups;

#define CUDA_CHECK(call)                                                    \
    do {                                                                    \
        cudaError_t err = call;                                             \
        if (err != cudaSuccess) {                                           \
            fprintf(stderr, "CUDA error %s:%d: %s\n", __FILE__, __LINE__,   \
                    cudaGetErrorString(err));                               \
            exit(1);                                                        \
        }                                                                   \
    } while (0)

__device__ double update_community_degree(int n_nodes, int *d_csr_adj_node, double *d_csr_weights, int *d_csr_node_offset, int *community, double *community_degree, int node_to_move, int old_community, int new_community) {
    for(int i = d_csr_node_offset[node_to_move]; i < d_csr_node_offset[node_to_move + 1]; i++){
        if(community[d_csr_adj_node[i]] == old_community){
            community_degree[old_community] -= d_csr_weights[i];
        }
        if(community[d_csr_adj_node[i]] == new_community){
            community_degree[new_community] += d_csr_weights[i];
        }
    }
    return 0.0;
}

__device__ double calculate_modularity_change(
        int n_nodes, 
        int m_edges, 
        double* d_total_weight, 
        int *d_csr_adj_node, 
        double *d_csr_weights, 
        int *d_csr_node_offset, 
        int *community, 
        double *community_degree, 
        int node_to_move, 
        int old_community, 
        int new_community) 
{
    double m = *d_total_weight;

    double delta_Q = 0.0;
    double k_i_in_old = 0.0, k_i_in_new = 0.0;
    double k_i = 0.0;
    double sum_tot_old = community_degree[old_community], sum_tot_new = community_degree[new_community];

    for(int i = d_csr_node_offset[node_to_move]; i < d_csr_node_offset[node_to_move + 1]; i++){
        k_i += d_csr_weights[i];
        if(community[d_csr_adj_node[i]] == old_community){
            k_i_in_old += d_csr_weights[i];       
        }
        if(community[d_csr_adj_node[i]] == new_community){
            k_i_in_new += d_csr_weights[i];    
        }
    }


    // change to printf for cuda debugging
    printf("Node %d: k_i_in_old = %f, k_i_in_new = %f, k_i = %f, sum_tot_old = %f, sum_tot_new = %f, m = %f\n", node_to_move, k_i_in_old, k_i_in_new, k_i, sum_tot_old, sum_tot_new, m);
    printf("%f This is internal changes\n", 2.0 * (k_i_in_new - k_i_in_old) / m);
    printf("sum_tot_old: %f sum_tot_new: %f k_i: %f\n", sum_tot_old, sum_tot_new, k_i);
    printf("%f This is degree changes\n", 2.0 * ((k_i * (sum_tot_old - sum_tot_new - k_i)) / (m * m)));

    delta_Q = 2.0 * (k_i_in_new - k_i_in_old) / m + 2.0 * ((k_i * (sum_tot_old - sum_tot_new - k_i)) / (m * m));
    return delta_Q;
}

__device__ int changed;

// kernel created for launch on only one block
// To allow multiple blocks, we need to use global memory to store changed flag
__global__ void louvain_kernel(
    int n_nodes,
    int m_edges,
    int* d_csr_ini_node,
    int *d_csr_adj_node,
    double *d_csr_weights,
    int *d_csr_node_offset,
    int *community,
    double *community_degree,
    int *d_vertex_locks,
    int *d_community_locks,
    double* d_total_weight ) 
{
    cg::grid_group grid = cg::this_grid();

    // Added /32 because only want one thread working per warp to avoid deadlock type situation
    int tid = blockIdx.x * (blockDim.x / 32) + threadIdx.x / 32;
    int nthreads = blockDim.x * gridDim.x / 32;
    if(threadIdx.x % 32 != 0) return;

    if(tid == 0) changed = 1;
    grid.sync();

    while(changed) {
        if(tid == 0) changed = 0;
        grid.sync();

        // Separate into sections and move
        // We discussed that vertex parallelism may be useful but came to a conclusion that sure it helps in 
        // the aspect of one lock but the other lock still needs to be taken.
        for(int i = (tid * m_edges) / nthreads; i < ((tid + 1) * m_edges) / nthreads; i++) {
            int node_to_move = d_csr_ini_node[i];
            int target_community_node = d_csr_adj_node[i];

            // take lock on nodes in order
            int first_node_lock = min(node_to_move, target_community_node);
            int second_node_lock = max(node_to_move, target_community_node);
            while (atomicCAS(&d_vertex_locks[first_node_lock], 0, 1) != 0) {}
            while (atomicCAS(&d_vertex_locks[second_node_lock], 0, 1) != 0) {}

            int initial_community = community[node_to_move];
            int target_community = community[target_community_node];

            if(initial_community != target_community) {
                // Take a lock on initial and target community in order of their ids to avoid deadlock
                int first_lock = min(initial_community, target_community);
                int second_lock = max(initial_community, target_community);

                // Wait until we can acquire both locks
                while (atomicCAS(&d_community_locks[first_lock], 0, 1) != 0) {}
                while (atomicCAS(&d_community_locks[second_lock], 0, 1) != 0) {}

                double delta_Q = calculate_modularity_change(n_nodes, m_edges, d_total_weight, d_csr_adj_node, d_csr_weights, d_csr_node_offset, community, community_degree, node_to_move, initial_community, target_community);
                printf("Node %d moving from community %d to %d gives delta Q = %f\n", node_to_move, initial_community, target_community, delta_Q);
                if(delta_Q > 1e-12) {
                    community[node_to_move] = target_community;
                    update_community_degree(n_nodes, d_csr_adj_node, d_csr_weights, d_csr_node_offset, community, community_degree, node_to_move, initial_community, target_community);
                    atomicExch(&changed, 1);
                }

                // Release the locks
                atomicExch(&d_community_locks[second_lock], 0);
                atomicExch(&d_community_locks[first_lock], 0);
            }

            // Release the node locks
            atomicExch(&d_vertex_locks[second_node_lock], 0);
            atomicExch(&d_vertex_locks[first_node_lock], 0);
        }

        grid.sync();
        // Print h_community for debugging
        if(tid == 0) {
            printf("Current community assignments:\n");
            for(int i = 0; i < n_nodes; i++) {
                printf("Node %d: Community %d\n", i, community[i]);
            }
        }

        grid.sync();
    }

}

__global__ void count_communities(
    int n_nodes,
    int *d_community, 
    int *d_community_present,
    int *d_num_communities)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int nthreads = blockDim.x * gridDim.x;
    
    for(int i = tid; i < n_nodes; i += nthreads) {
        int comm = d_community[i];
        atomicExch(&d_community_present[comm], 1);
    }

    __shared__ int community_count_in_block;
    if(threadIdx.x == 0) {
        community_count_in_block = 0;
    }
    __syncthreads();

    for(int i = tid; i < n_nodes; i += nthreads) {
        if(d_community_present[i] == 1) {
            atomicAdd(&community_count_in_block, 1);
        }
    }

    if(threadIdx.x == 0) atomicAdd(d_num_communities, community_count_in_block);    
}

void louvain_cuda(vector<vector<pair<int, double>>> &adj, int n_nodes, int m_edges) {
    int n = adj.size();

    // Create CSR representation
    int* h_csr_ini_node = new int[m_edges];
    int* h_csr_adj_node = new int[m_edges];
    double* h_csr_weights = new double[m_edges];
    int* h_csr_node_offset = new int[n_nodes + 1];
    int* h_community = new int[n_nodes];
    double* h_total_weight = new double; *h_total_weight = 0.0; // actually double total weight
    double* h_community_degree = new double[n_nodes];
    for (int i = 0; i < n_nodes; i++) {
        h_community_degree[i] = 0.0;
    }

    int edge_counter = 0;
    for (int i = 0; i < n_nodes; i++) {
        h_csr_node_offset[i] = edge_counter;
        h_community[i] = i;
        for (const auto &edge : adj[i]) {
            h_csr_ini_node[edge_counter] = i;
            h_csr_adj_node[edge_counter] = edge.first;
            h_csr_weights[edge_counter] = edge.second;
            edge_counter++;
            h_community_degree[i] += edge.second;
            h_total_weight[0] += edge.second;
        }
    }
    h_csr_node_offset[n_nodes] = m_edges;

    // Allocate device memory
    int *d_csr_ini_node, *d_csr_adj_node, *d_csr_node_offset, *d_community;
    int *d_community_locks, *d_vertex_locks;
    double *d_total_weight, *d_csr_weights, *d_community_degree;
    
    cudaMalloc(&d_csr_ini_node, m_edges * sizeof(int));
    cudaMalloc(&d_csr_adj_node, m_edges * sizeof(int));
    cudaMalloc(&d_csr_weights, m_edges * sizeof(double));
    cudaMalloc(&d_csr_node_offset, (n_nodes + 1) * sizeof(int));
    cudaMalloc(&d_community, n_nodes * sizeof(int));
    cudaMalloc(&d_community_degree, n_nodes * sizeof(double));
    cudaMalloc(&d_vertex_locks, n_nodes * sizeof(int));
    cudaMemset(d_vertex_locks, 0, n_nodes * sizeof(int));
    cudaMalloc(&d_community_locks, n_nodes * sizeof(int));
    cudaMemset(d_community_locks, 0, n_nodes * sizeof(int));
    cudaMalloc(&d_total_weight, sizeof(double));

    cudaMemcpy(d_csr_ini_node, h_csr_ini_node, m_edges * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_csr_adj_node, h_csr_adj_node, m_edges * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_csr_weights, h_csr_weights, m_edges * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_csr_node_offset, h_csr_node_offset, (n_nodes + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_community, h_community, n_nodes * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_community_degree, h_community_degree, n_nodes * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_total_weight, h_total_weight, sizeof(double), cudaMemcpyHostToDevice);

    int* d_community_present, *d_num_communities;
    cudaMalloc(&d_community_present, n_nodes * sizeof(int));
    cudaMemset(d_community_present, 0, n_nodes * sizeof(int));
    cudaMalloc(&d_num_communities, sizeof(int));
    cudaMemset(d_num_communities, 0, sizeof(int));

    // Cooperative Launch
    void* kernelArgs[] = {
        (void*)&n_nodes,
        (void*)&m_edges,
        (void*)&d_csr_ini_node,
        (void*)&d_csr_adj_node,
        (void*)&d_csr_weights,
        (void*)&d_csr_node_offset,
        (void*)&d_community,
        (void*)&d_community_degree,
        (void*)&d_vertex_locks,
        (void*)&d_community_locks,
        (void*)&d_total_weight
    };
    cudaLaunchCooperativeKernel((void*)louvain_kernel, 32, 512, kernelArgs);

    // Non cooperative launch
    // louvain_kernel<<<1024, 1>>>(
    //     n_nodes, 
    //     m_edges,
    //     d_csr_ini_node, 
    //     d_csr_adj_node, 
    //     d_csr_weights, 
    //     d_csr_node_offset, 
    //     d_community, 
    //     d_community_degree,
    //     d_vertex_locks,  
    //     d_community_locks, 
    //     d_total_weight);
    cudaDeviceSynchronize();
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy data from device to host
    cudaMemcpy(h_csr_ini_node, d_csr_ini_node, m_edges * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_csr_adj_node, d_csr_adj_node, m_edges * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_csr_weights, d_csr_weights, m_edges * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_csr_node_offset, d_csr_node_offset, (n_nodes + 1) * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_community, d_community, n_nodes * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_community_degree, d_community_degree, n_nodes * sizeof(double), cudaMemcpyDeviceToHost);

    count_communities<<<32, 512>>>(n_nodes, d_community, d_community_present, d_num_communities);
    cudaDeviceSynchronize();
    int h_num_communities;
    cudaMemcpy(&h_num_communities, d_num_communities, sizeof(int), cudaMemcpyDeviceToHost);
    cout << "Number of communities: " << h_num_communities << endl;

    thrust::device_ptr<int> dev_ptr = thrust::device_pointer_cast(d_community_present);
    thrust::inclusive_scan(dev_ptr, dev_ptr + n_nodes, dev_ptr);
    for(int i = 0; i < n_nodes; i++) {
        int temp;
        cudaMemcpy(&temp, d_community_present + i, sizeof(int), cudaMemcpyDeviceToHost);
        cout << "Community presence at (Prefix Summed)" << i << " is " << temp << endl;
    }

    // Output the communities
    for (int i = 0; i < n_nodes; i++) {
        cout << "Partition " << i << " is in community " << h_community[i] << endl;
    }

    // Change edges based on the new community indexes
    // Change the community array based on the new indexes
    // May want to put edges that become between the same two partitions

    n_nodes = h_num_communities;
    // aggregate_graph<<<32, 512>>>(d_csr_ini_node,
    //     d_csr_adj_node,
    //     d_csr_weights,
    //     d_csr_node_offset,
    //     d_community,
    //     n_nodes,
    //     m_edges);

    // Free host memory
    delete[] h_csr_adj_node;
    delete[] h_csr_weights;
    delete[] h_csr_node_offset;
    delete[] h_community;
    delete[] h_community_degree;
}

int main() {
    int n_nodes, m_edges;
    cin >> n_nodes >> m_edges;
    vector<vector<pair<int, double>>> adj(n_nodes);
    for (int i = 0; i < m_edges; i++) {
        int u, v;
        double w;
        cin >> u >> v >> w;
        // Imagine using emplace back instead of push back
        adj[u].emplace_back(v, w);
        adj[v].emplace_back(u, w);
    }

    louvain_cuda(adj, n_nodes, m_edges * 2);
    return 0;
}