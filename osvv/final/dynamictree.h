#include "stdbool.h"

typedef struct dynamic_tree_node {
    int s, my_s, on, id;
    bool flip, my_flip;
    struct dynamic_tree_node* l, r, p;

} dynamic_node_t;


void init(dynamic_node_t* node, int c, int i);

bool isroot(dynamic_node_t* node);

/* If this node is flipped, we unflip it, and push the change
   down the tree, so that it represents the same thing. */
void normalize(dynamic_node_t* node);

/* The tree structure has changed in the vicinity of this node
   (for example, if this node is linked to a different left
   child in a rotation).  This function fixes up the data fields
   in the node to maintain invariants. */
void update(dynamic_node_t* node);



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
int toggle(dynamic_node_t* p);

/* this returns the id of the node that is the root of the tree containing p */
int rootid(dynamic_node_t* p);
