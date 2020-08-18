#include "stdbool.h"
#include "types.h"

typedef struct dynamic_tree_node {
    long s, my_s, on, id;
    bool flip, my_flip;
    long delmin, delcost;
    arc* edge;
    struct dynamic_tree_node* l, r, p;

} dynamic_node_t;

dynamic_node_t* nodes;


typedef struct dynamic_tree_path {
    dynamic_node_t* root;
    dynamic_node_t* head, tail;

} dynamic_path_t;


void init(dynamic_node_t* node, long c, long i);

bool isroot(dynamic_node_t* node);

/* If this node is flipped, we unflip it, and push the change
   down the tree, so that it represents the same thing. */
void normalize(dynamic_node_t* node);

/* The tree structure has changed in the vicinity of this node
   (for example, if this node is linked to a different left
   child in a rotation).  This function fixes up the data fields
   in the node to maintain invariants. */
void update(dynamic_node_t* node;



void rotR (dynamic_node_t* p);

void rotL (dynamic_node_t* p);

void splay(dynamic_node_t* p);

/* This makes node q the root of the virtual tree, and also q is the
   leftmost node in its splay tree */
void expose(dynamic_node_t* q);

/* assuming p and q are nodes in different trees and
   that p is a root of its tree, this links p to q */
void link(dynamic_node_t* p, dynamic_node_t* q);

    /* Toggle all the edges on the path from p to the root
       return the count after - count before */
// long toggle(dynamic_node_t* p);

/* this returns the id of the node that is the root of the tree containing p */
long rootid(dynamic_node_t* p);



/* each node in dynamic_path_t stores the cost of the edge to its parent,
   except for the leftmost node, i.e., the head of the path*/
// dt_path_t* concatenate(dynamic_path_t* p, dynamic_path_t* q, arc* edge);
dynamic_node_t* concatenate(dynamic_node_t* p, dynamic_node_t* q, dynamic_node_t* r);
dynamic_node_t* before(dynamic_node_t* v);
dynamic_node_t* after(dynamic_node_t* v);

long pMinCost(dynamic_path_t* p);
long nMinCost(dynamic_node_t* v);
long nCost(dynamic_node_t* v);
void pUpdate(dynamic_node_t* p, long x);
void cut(dynamic_node_t* v);
void findPath(dynamic_path_t* p, long* a, long* b, long* cost);
