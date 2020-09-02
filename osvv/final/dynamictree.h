#include "stdbool.h"
#include "types.h"

#ifndef FILE_DYNAMICTREE_SEEN
#define FILE_DYANMICTREE_SEEN

typedef struct dynamic_tree_node {
    long id;
    long delmin, delcost;
    arc* edge; // edge going out from this vertex
    struct dynamic_tree_node *left, *right, *parent;
} dynamic_node_t;


typedef struct dynamic_tree {
    long sz;
    dynamic_node_t* nodes;
    long cur_node;

} dynamic_tree_t;



dynamic_tree_t* init(long nodes, long start_node);

void initNode(dynamic_node_t* node, long i);

void cleanUp(dynamic_tree_t* dTree);

bool isroot(dynamic_node_t* node);



void rotR (dynamic_node_t* p);

void rotL (dynamic_node_t* p);

void splay(dynamic_node_t* p);

/* This makes node q the root of the virtual tree, and also q is the
   leftmost node in its splay tree */
void expose(dynamic_node_t* q);

/* assuming p and q are nodes in different trees and
   that p is a root of its tree, this links p to q */
void link(dynamic_tree_t* dTree, long pid, long qid, arc* edge);


/* this returns the id of the node that is the root of the tree containing p */
long root(dynamic_tree_t* dTree, long pid);




long before(dynamic_tree_t* dTree, long vid);
long after(dynamic_tree_t* dTree, long vid);

// long pMinCost(dynamic_path_t* p);
long nMinCost(dynamic_tree_t* dTree, long vid);
long nCost(dynamic_tree_t* dTree, long vid);
void pUpdate(dynamic_tree_t* dTree, long pid, long x);
void cut(dynamic_tree_t* dTree, long vid);
void cutEdge(dynamic_tree_t* dTree, long vid);
void findPath(dynamic_tree_t* dTree, long pid, long* a, long* b, long* cost);


#endif
