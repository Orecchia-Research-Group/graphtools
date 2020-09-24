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
    dynamic_node_t* d_nodes;
    node* nodes;
    dynamic_node_t* d_cur_node;  // the root of the tree containing source
    node* cur_node;
} dynamic_tree_t;



dynamic_tree_t* dec_init(long num_nodes, node* nodes, node* start_node);

void initNode(dynamic_node_t* node, long i);

void cleanUp(dynamic_tree_t* dTree);

dynamic_node_t* to_d_node(dynamic_tree_t* dTree, node* p);

node* to_node(dynamic_tree_t* dTree, dynamic_node_t* p);

bool isroot(dynamic_node_t* node);



void rotR (dynamic_node_t* p);

void rotL (dynamic_node_t* p);

void splay(dynamic_node_t* p);

/* This makes node q the root of the virtual tree, and also q is the
   leftmost node in its splay tree */
void expose(dynamic_node_t* q);

/* assuming p and q are nodes in different trees and
   that p is a root of its tree, this links p to q */
void link(dynamic_tree_t* dTree, dynamic_node_t* p, dynamic_node_t* q, arc* edge);


/* this returns the id of the node that is the root of the tree containing p */
dynamic_node_t* d_root(dynamic_node_t* p);
node* root(dynamic_tree_t* dTree, node* p);




dynamic_node_t* d_before(dynamic_node_t* p);
dynamic_node_t* d_after(dynamic_node_t* p);
node* before(dynamic_tree_t* dTree, node* p);
node* after(dynamic_tree_t* dTree, node* p);

long nMinCost(dynamic_node_t* p);
long nCost(dynamic_node_t* p);
void pUpdate(dynamic_node_t* p, long x);

void d_cut(dynamic_tree_t* dTree, dynamic_node_t* p);
void cut(dynamic_tree_t* dTree, node* p);

void d_cutEdge(dynamic_tree_t* dTree, dynamic_node_t* p);
void cutEdge(dynamic_tree_t* dTree, node* p);
void findPath(dynamic_tree_t* dTree, node* p, node** a, node** b, long* cost);


#endif
