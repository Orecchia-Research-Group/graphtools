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















void rotR (dynamic_node_t* p) {
    dynamic_node_t* q = p->p;
    dynamic_node_t* r = q->p;
    normalize(q);
    normalize(p);
    if ((q->l=p->r) != NULL) {
        q->l->p = q;
    }
    p->r = q;
    q->p = p;
    if ((p->p=r) != NULL) {
        if (r->l == q) r->l = p;
        else if (r->r == q) r->r = p;
    }
    update(q);
}

void rotL (dynamic_node_t* p) {
    dynamic_node_t* q = p->p;
    dynamic_node_t* r = q->p;
    normalize(p);
    normalize(q);
    if ((q->r=p->l) != NULL) {
        q->r->p = q;
    }
    p->l = q;
    q->p = p;
    if ((p->p=r) != NULL) {
        if (r->l == q) r->l = p;
        else if (r->r == q) r->r = p;
    }
    update(q);
}

void splay(dynamic_node_t* p) {
    while (!isroot(p)) {
        dynamic_node_t* q = p->p;
        if (isroot(q)) {
            if (q->l == p)
                rotR(p);
            else
                rotL(p);
        }
        else {
            dynamic_node_t* r = q->p;
            if (r->l == q) {
                if (q->l == p) {
                    rotR(q);
                    rotR(p);
                }
                else {
                    rotL(p);
                    rotR(p);
                }
            }
            else {
                if (q->r == p) {
                    rotL(q);
                    rotL(p);
                }
                else {
                    rotR(p);
                    rotL(p);
                }
            }
        }
    }
    normalize(p); // only useful if p was already a root.
    update(p);    // only useful if p was not already a root
}

/* This makes node q the root of the virtual tree, and also q is the
   leftmost node in its splay tree */
void expose(dynamic_node_t* q) {
    dynamic_node_t* r = NULL;
    for (dynamic_node_t* p=q; p != NULL; p=p->p) {
        splay(p);
        p->l = r;
        update(p);
        r = p;
    }
    splay(q);
}

/* assuming p and q are nodes in different trees and
   that p is a root of its tree, this links p to q */
void link(dynamic_node_t* p, dynamic_node_t* q) {
    expose(p);
    if (p->r != NULL) {
        // p is not a root. Error
        return;
    }
    p->p = q;
}

    /* Toggle all the edges on the path from p to the root
       return the count after - count before */
int toggle(dynamic_node_t* p) {
    expose(p);
    int before = p->on;
    p->flip = !p->flip;
    normalize(p);
    int after = p->on;
    return after - before;
}

/* this returns the id of the node that is the root of the tree containing p */
int rootid(dynamic_node_t* p) {
    expose(p);
    while(p->r != NULL) p = p->r;
    splay(p);
    return p->id;
}
