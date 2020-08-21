#include "stdbool.h"
#include "types.h"

#ifndef FILE_DYNAMICTREE_SEEN
#define FILE_DYANMICTREE_SEEN

typedef struct dynamic_tree_node {
    long id;
    long delmin, delcost;
    arc* edge; // edge going out from this vertex
    struct dynamic_tree_node *l, *r, *p;
} dynamic_node_t;


typedef struct dynamic_tree {
    long sz;
    dynamic_node_t* nodes;

} dynamic_tree_t;

dynamic_tree_t dTree;



void init(long nodes);

void initNode(dynamic_node_t* node, long i);

void cleanUp();

bool isroot(dynamic_node_t* node);



void rotR (dynamic_node_t* p);

void rotL (dynamic_node_t* p);

void splay(dynamic_node_t* p);

/* This makes node q the root of the virtual tree, and also q is the
   leftmost node in its splay tree */
void expose(dynamic_node_t* q);

/* assuming p and q are nodes in different trees and
   that p is a root of its tree, this links p to q */
void link(long pid, long qid, arc* edge);

    /* Toggle all the edges on the path from p to the root
       return the count after - count before */
// long toggle(dynamic_node_t* p);

/* this returns the id of the node that is the root of the tree containing p */
long root(long pid);




long before(long vid);
long after(long vid);

// long pMinCost(dynamic_path_t* p);
long nMinCost(long vid, dynamic_node_t** rootptr);
long nCost(long vid);
void pUpdate(long pid, long x);
void cut(long vid);
void cutEdge(long vid);
void findPath(long pidi, long* a, long* b, long* cost);


#endif
