The modularity does not change when the graph is aggregated.
My implementation in dynamic_louvain_dynamic_frontier involves first doing the static louvain algorithm.
After this I deaggregate the graph. And apply the dynamic_frontier algorithm.
In this case I consider each partition as one node. A partition is simply a group of nodes moving together as one group.

Optimizations to do in the Parallel Setting:
- Modularity recomputation is not necessary
- Go through harshitha's code
- Need to parallelly move nodes even in louvain
- Need to represent the graph in a format that is easily changeable

- When you move u to v, communities not directly connected are not affected?

This README is DEPRECATED by me :)
All notes being taken on personal google docs