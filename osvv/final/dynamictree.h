#include "stdbool.h"
#include "types.h"

#ifndef FILE_DYNAMICTREE_SEEN
#define FILE_DYANMICTREE_SEEN

//#define DEBUG

typedef struct dynamic_tree_node {
    long id;
    long delmin, delcost;
    arc* edge; // edge going out from this vertex
    struct dynamic_tree_node *left, *right, *parent;

#ifdef DEBUG
    long rot_cnt;
    long splay_cnt;
    long expose_cnt;
#endif
} dynamic_node_t;


typedef struct dynamic_tree {
    long sz;
    dynamic_node_t* d_nodes;
    node* nodes;
    dynamic_node_t* d_cur_node;  // the root of the tree containing source
    node* cur_node;
    dynamic_node_t* d_source;
    node* source;

#ifdef DEBUG
    long link_cnt;
    long cut_cnt;
#endif
} dynamic_tree_t;

void dt_print_op_stat(dynamic_tree_t* dTree);

dynamic_tree_t* dt_init(long num_nodes, node* nodes, node* start_node);

void dt_initNode(dynamic_node_t* node, long i);

void dt_cleanUp(dynamic_tree_t* dTree);

dynamic_node_t* dt_to_d_node(dynamic_tree_t* dTree, node* p);

node* dt_to_node(dynamic_tree_t* dTree, dynamic_node_t* p);

bool dt_isroot(dynamic_node_t* node);



void dt_rotR (dynamic_node_t* p);

void dt_rotL (dynamic_node_t* p);

void dt_splay(dynamic_node_t* p);

/* This makes node q the root of the virtual tree, and also q is the
   leftmost node in its splay tree */
void dt_expose(dynamic_node_t* q);

/* assuming p and q are nodes in different trees and
   that p is a root of its tree, this links p to q,
   and returns true. If q is the tree rooted at p,
   this removes a cycle and returns false*/
bool dt_d_link(dynamic_tree_t* dTree, dynamic_node_t* p, dynamic_node_t* q, arc* edge);
bool dt_link(dynamic_tree_t* dTree, node* p, node* q, arc* edge);


/* this returns the id of the node that is the root of the tree containing p */
dynamic_node_t* dt_d_root(dynamic_node_t* p);
node* dt_root(dynamic_tree_t* dTree, node* p);




dynamic_node_t* dt_d_before(dynamic_node_t* p);
dynamic_node_t* dt_d_after(dynamic_node_t* p);
node* dt_before(dynamic_tree_t* dTree, node* p);
node* dt_after(dynamic_tree_t* dTree, node* p);

long dt_nMinCost(dynamic_node_t* p);
long dt_nCost(dynamic_node_t* p);
void dt_pUpdate(dynamic_node_t* p, long x);

void dt_d_cut(dynamic_tree_t* dTree, dynamic_node_t* p);
void dt_cut(dynamic_tree_t* dTree, node* p);

void dt_d_cutEdge(dynamic_tree_t* dTree, dynamic_node_t* p);
void dt_cutEdge(dynamic_tree_t* dTree, node* p);
void dt_findPath(dynamic_tree_t* dTree, node** a, node** b, long* cost);

long dt_dfs(dynamic_tree_t* p);


#endif
