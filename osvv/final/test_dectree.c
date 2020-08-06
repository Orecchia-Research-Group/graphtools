// #include "types.h"
#include "dectree.h"
#include "stdlib.h"
#include "stddef.h"
#include "stdbool.h"
#include "stdio.h"
#include "assert.h"


int main() {


    decInit(4,5);
    assert(path(0) != NULL);
    assert(path(1) != NULL);
    // nodes are 0,1,2,3
    // edges are (0, 1, 1), (0, 2, 3), (2, 1, 2), (2, 3, 1), (3, 1, 1)
    arc* edges = calloc(5, sizeof(arc));
    edges[0].resCap = 1;
    edges[1].resCap = 3;
    edges[2].resCap = 2;
    edges[3].resCap = 1;
    edges[4].resCap = 1;



    long* match_a = calloc(3, sizeof(long));
    long* match_b = calloc(3, sizeof(long));
    long* match_cost = calloc(3, sizeof(long));







    dt_path_t* p0 = path(0);
    dt_path_t* p1 = path(1);
    assert(p0 != NULL && p1 != NULL);
    p0 = concatenate(p0, p1, &edges[0]);
    long cur = tail(p0);
    long curh = head(p0);
    printf("%ld %ld\n", curh, cur);
    assert(cur == 1);

    findPath(p0, &match_a[0], &match_b[0], &match_cost[0]);
    printf("match1 = %ld %ld %ld\n", match_a[0], match_b[0], match_cost[0]);

    cur = tail(p0);
    assert(cur == 0);
    dt_path_t* p2 = path(2);
    p0 = concatenate(p0, p2, &edges[1]);
    p1 = path(1);
    p0 = concatenate(p0, p1, &edges[2]);
    cur = tail(p0);
    assert(cur == 1);

    findPath(p0, &match_a[1], &match_b[1], &match_cost[1]);

    cur = tail(p0);
    assert(cur == 2);
    dt_path_t* p3 = path(3);
    p0 = concatenate(p0, p3, &edges[3]);
    p1 = path(1);
    p0 = concatenate(p0, p1, &edges[4]);
    cur = tail(p0);
    assert(cur == 1);

    findPath(p0, &match_a[2], &match_b[2], &match_cost[2]);

    cur = tail(p0);
    assert(cur == 0);


    return 0;
}
