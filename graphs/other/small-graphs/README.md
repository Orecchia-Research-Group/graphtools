#hiprtest graphs
The graph testhipr.eg2 and testhipr2.eg2 test the cutfind function.

The following characteristics are tested for:
	1. If the algoritm goes along reverse edges, resulting in an infinite loop. Tested using a cycle of length 3 and an edge that is only incident to the sink.
	2. If the algorithm goes along 0 flow edges. Tested using multiple edges in the same paths which will have resCap of 0 eventually.

A good method to check if the edge capacities are being update appropriately: comparing the flow outgoing from the source and flow ingoing into the sink, which should be equal. 
 
