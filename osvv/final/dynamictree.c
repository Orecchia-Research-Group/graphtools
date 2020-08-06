#include "dynamictree.h"




void init(dynamic_node_t* node, int c, int i) {
    node->id = i;
    node->s = node->my_s = c;
    node->on = 0;
    node->l = node->r = node->p = NULL;
    node->flip = node->my_flip = false;
}

bool isroot(dynamic_node_t* node) {
       return (node->p==NULL) || (node->p->l != node && node->p->r != node);
}

/* If this node is flipped, we unflip it, and push the change
   down the tree, so that it represents the same thing. */
void normalize(dynamic_node_t* node) {
    if (node->flip) {
        node->flip = false;
        node->on = node->s - node->on;
        node->my_flip = !node->my_flip;
        if (node->l != NULL) node->l->flip = !node->l->flip;
        if (node->r != NULL) node->r->flip = !node->r->flip;
    }
}

/* The tree structure has changed in the vicinity of this node
   (for example, if this node is linked to a different left
   child in a rotation).  This function fixes up the data fields
   in the node to maintain invariants. */
void update(dynamic_node_t* node) {
    node->s = node->my_s;
    node->on = (node->my_flip)?node->my_s:0;
    if (node->l != NULL) {
        node->s += node->l->s;
        if (node->l->flip)
            node->on += node->l->s - node->l->on;
        else
            node->on += node->l->on;
    }
    if (node->r != NULL) {
        node->s += node->r->s;
        if (node->r->flip)
            node->on += node->r->s - node->r->on;
        else node->on += node->r->on;
    }
}
