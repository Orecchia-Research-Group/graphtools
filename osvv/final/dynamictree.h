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
