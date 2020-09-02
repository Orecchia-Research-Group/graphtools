#include "dynamictree.h"
#include "stdlib.h"
#include "stddef.h"
#include "stdbool.h"
#include "stdio.h"
#include "assert.h"


int main() {
    // nodes are 0,1,2,3,4
    // edges are (0, 1, 1), (0, 2, 3), (2, 1, 2), (2, 3, 1), (3, 1, 1), (1, 4, 4)
    arc* edges = calloc(6, sizeof(arc));
    edges[0].cap = 1;
    edges[1].cap = 3;
    edges[2].cap = 2;
    edges[3].cap = 1;
    edges[4].cap = 1;
    edges[5].cap = 4;
    edges[0].resCap = 0;
    edges[1].resCap = 0;
    edges[2].resCap = 0;
    edges[3].resCap = 0;
    edges[4].resCap = 0;
    edges[5].resCap = 0;

    dynamic_tree_t* dTree = init(5, 0);



    long* match_a = calloc(3, sizeof(long));
    long* match_b = calloc(3, sizeof(long));
    long* match_cost = calloc(3, sizeof(long));

    link(dTree, 0, 1, &edges[0]);
    printf("link done\n");
    assert(dTree->cur_node == 1);
    link(dTree, 1, 4, &edges[5]);
    printf("link done\n");
    assert(dTree->cur_node == 4);
    findPath(dTree, 0, &match_a[0], &match_b[0], &match_cost[0]);
    printf("%ld %ld %ld\n", match_a[0], match_b[0], match_cost[0]);
    assert(dTree->cur_node == 0);

    link(dTree, 0, 2, &edges[1]);
    printf("link done\n");
    assert(dTree->cur_node == 2);
    link(dTree, 2, 1, &edges[2]);
    printf("link done\n");
    assert(dTree->cur_node == 4);
    findPath(dTree, 0, &match_a[1], &match_b[1], &match_cost[1]);
    printf("%ld %ld %ld\n", match_a[1], match_b[1], match_cost[1]);
    assert(dTree->cur_node == 2);

    link(dTree, 2, 3, &edges[3]);
    printf("link done\n");
    assert(dTree->cur_node == 3);
    link(dTree, 3, 1, &edges[4]);
    printf("link done\n");
    assert(dTree->cur_node == 4);
    findPath(dTree, 0, &match_a[2], &match_b[2], &match_cost[2]);
    printf("%ld %ld %ld\n", match_a[2], match_b[2], match_cost[2]);
    assert(dTree->cur_node == 0);

    cleanUp(dTree);


    return 0;
}
