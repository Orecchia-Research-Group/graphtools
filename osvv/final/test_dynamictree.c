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
    edges[0].resCap = 1;
    edges[1].resCap = 3;
    edges[2].resCap = 2;
    edges[3].resCap = 1;
    edges[4].resCap = 1;
    edges[5].resCap = 4;

    init(5);



    long* match_a = calloc(3, sizeof(long));
    long* match_b = calloc(3, sizeof(long));
    long* match_cost = calloc(3, sizeof(long));

    link(0, 1, &edges[0]);
    printf("link done\n");
    link(1, 4, &edges[5]);
    printf("link done\n");
    findPath(0, &match_a[0], &match_b[0], &match_cost[0]);
    printf("%ld %ld %ld\n", match_a[0], match_b[0], match_cost[0]);

    link(0, 2, &edges[1]);
    printf("link done\n");
    link(2, 1, &edges[2]);
    printf("link done\n");
    findPath(0, &match_a[1], &match_b[1], &match_cost[1]);
    printf("%ld %ld %ld\n", match_a[1], match_b[1], match_cost[1]);

    link(2, 3, &edges[3]);
    printf("link done\n");
    link(3, 1, &edges[4]);
    printf("link done\n");
    findPath(0, &match_a[2], &match_b[2], &match_cost[2]);
    printf("%ld %ld %ld\n", match_a[2], match_b[2], match_cost[2]);



    // link(&nodes[5], &nodes[0]);
    // link(&nodes[0], &nodes[6]);
    // findPath(&nodes[5], &match_a[0], &match_b[0], &match_cost[0]);
    //
    //
    // assert(root(&nodes[5]) == &nodes[5]);
    // link(&nodes[5], &nodes[1]);
    // link(&nodes[1], &nodes[2]);
    // link(&nodes[2], &nodes[6]);
    // findPath(&nodes[5], &match_a[1], &match_b[1], &match_cost[1]);
    //
    // assert(root(&nodes[5]) == &nodes[1]);
    // link(&nodes[1], &nodes[3]);
    // link(&nodes[3], &nodes[4]);
    // link(&nodes[4], &nodes[6]);
    // findPath(&nodes[5], &match_a[2], &match_b[2], &match_cost[2]);

    return 0;
}
