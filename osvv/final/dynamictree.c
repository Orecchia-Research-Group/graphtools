#include "dynamictree.h"

#define min(a,b) ((a) < (b)) ? (a) : (b)
#define max(a,b) ((a) > (b)) ? (a) : (b)


void init(dynamic_node_t* node, long c, long i) {
    node->id = i;
    node->s = node->my_s = c;
    node->on = 0;
    node->l = node->r = node->p = NULL;
    node->flip = node->my_flip = false;
    node->delmin = node->delcost = 0;
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
    // normalize(q);
    // normalize(p);

    // update cost
    long p_delcost = p->delcost + q->delcost;
    long q_delcost = -p->delcost;
    if (p->r != NULL) {
        p->r->delcost = p->delcost + p->r->delcost;
    }
    long q_delmin = max(0, p->r->delmin - pr_delcost);
    if (q->r != NULL) {
        q_delmin = max(q_delmin, q->r->delmin - q->r->delcost);
    }
    long p_delmin = max(0, q_delmin - q_delcost);
    if (p->l != NULL) {
        p_delmin = max(p_delmin, p->l->delmin - p->l->delcost);
    }


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
    // normalize(p);
    // normalize(q);

    // update cost
    long p_delcost = p->delcost + q->delcost;
    long q_delcost = -p->delcost;
    if (p->l != NULL) {
        p->l->delcost = p->delcost + p->l->delcost;
    }
    long q_delmin = max(0, p->r->delmin - pr_delcost);
    if (q->l != NULL) {
        q_delmin = max(q_delmin, q->l->delmin - q->l->delcost);
    }
    long p_delmin = max(0, q_delmin - q_delcost);
    if (p->r != NULL) {
        p_delmin = max(p_delmin, p->r->delmin - p->r->delcost);
    }


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

        dynamic_node_t* s = p->l;
        if (s != NULL) {
            s->delcost += p->delcost;
        }

        p->l = r;


        p->delmin = 0;
        if (r != NULL) {
            r->delcost -= p->delcost;
            p->delmin = max(p->delmin, r->delmin - r->delcost);
        }

        if (p->r != NULL) {
            p->delmin = max(p->delmin, p->r->delmin - p->r->delcost);
        }


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
// long toggle(dynamic_node_t* p) {
//     expose(p);
//     long before = p->on;
//     p->flip = !p->flip;
//     normalize(p);
//     long after = p->on;
//     return after - before;
// }

/* this returns the id of the node that is the root of the tree containing p */
long rootid(dynamic_node_t* p) {
    expose(p);
    while(p->r != NULL) p = p->r;
    splay(p);
    return p->id;
}








// cost implementation

// dt_path_t* concatenate(dynamic_path_t* p, dynamic_path_t* q, arc* edge) {
//     dynamic_node_t* phead = p->head, ptail = p->tail, qhead = q->head, qtail = q->tail;
//     dynamic_node_t* proot = p->root, qroot = q->root;
//     // link(proot, qhead, edge);
//
//     // update vals for ptail
//     long min_cost = (ptail->l != NULL) ? minCost(ptail->l) : 0x3f3f3f3f3f3f3f3f;
//     ptail->delcost += edge->resCap;
//     long cost = nCost(ptail);
//     ptail->delmin = cost - min(cost, min_cost);
//
//
//     // structural change
//     proot->p = qhead;
//     qhead->l = proot;
//     p->root = qroot;
//     ptail->edge = edge;
//     p->tail = qtail;
//     free(q);
//
//     // update delmin
//     for (dynamic_node_t* u = ptail; u != NULL; u = u->p) {
//
//     }
// }

dynamic_node_t* concatenate(dynamic_node_t* p, dynamic_node_t* q, dynamic_node_t* r) {

}


dynamic_node_t* before(dynamic_node_t* v) {
    if (v->l != NULL) {
        dynamic_node_t* u = v->l;
        for (;u->r;) {
            u = u->r;
        }
        return u;
    }
    else if (v->p != NULL && v->p->r == v) {
        return v->p;
    }
    else {
        return NULL; // v is the head of the path, the leftmost node of the splay tree
    }
}



dynamic_node_t* after(dynamic_node_t* v) {
    if (v->r != NULL) {
        dynamic_node_t* u = v->r;
        for (;u->l;) {
            u = u->l;
        }
        return u;
    }
    else if (v->p != NULL && v->p->l == v) {
        return v->p;
    }
    else {
        return NULL; // v is the tail of the path, the rightmost node of the splay tree
    }
}


long pMinCost(dynamic_path_t* p) {
    return nCost(p->root) - p->root->delmin;
}


long nMinCost(dynamic_node_t* v) {
    // TODO: splay
    return nCost(v) - v->delmin;
}


long nCost(dynamic_node_t* v) {
    // TODO: splay
    long cost = 0;
    for(dynamic_node_t* u = v; u != NULL; u = u->p) {
        cost += u->delcost;
    }
    return cost;
}


void pUpdate(dynamic_node_t* p, long x) {
    expose(p);
    p->delcost += x;
    // dynamic_node_t* proot = p->root;
    // proot->delcost += x;
    // long cost = nCost(proot);
    // long min_cost = 0x3f3f3f3f3f3f3f3f;
    // if (proot->l != NULL) {
    //     min_cost = min(min_cost, nMinCost(proot->l));
    // }
    // if (proot->r != NULL) {
    //     min_cost = min(min_cost, nMinCost(proot->r));
    // }
    // proot->delmin = cost - min(cost, min_cost);
}


void cut(dynamic_node_t* v) {
    splay(v);
    v->l->p = NULL;
    v->r->p = NULL;
    v->l = v->r = NULL;

    // update cost
    v->l->delcost += v->delcost;
    v->r->delcost += v->delcost;
}


void findPath(dynamic_path_t* p, long* a, long* b, long* cost) {


}
