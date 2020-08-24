#include "dynamictree.h"
#include "stdlib.h"
#include "stddef.h"
#include "stdio.h"

#define min(a,b) ((a) < (b)) ? (a) : (b)
#define max(a,b) ((a) > (b)) ? (a) : (b)

#define inf (0x3f3f3f3f3f3f3f3fL)


void init(long nodes, long start_node) {
    dTree.sz = nodes;
    dTree.nodes = calloc(nodes, sizeof(dynamic_node_t));
    for (long i = 0; i < dTree.sz; i++) { // (TODO: int or long)
        initNode(&dTree.nodes[i], i);
    }
    dTree.cur_node = start_node;
}

void initNode(dynamic_node_t* node, long i) {
    node->id = i;
    node->l = node->r = node->p = NULL;
    node->delmin = node->delcost = 0;
    node->edge = NULL;
}

void cleanUp() {

}

bool isroot(dynamic_node_t* node) {
    return (node->p==NULL) || (node->p->l != node && node->p->r != node);
}





/*             q              p
 *           /  \           /  \
 *         p     c   ->   a     q
 *       /  \                 /  \
 *     a     b              b     c
 */
void rotR (dynamic_node_t* p) {
    dynamic_node_t* q = p->p;
    dynamic_node_t* r = q->p;
    dynamic_node_t* a = p->l;
    dynamic_node_t* b = p->r;
    dynamic_node_t* c = q->r;
    long p_delcost = p->delcost + q->delcost;
    long q_delcost = -p->delcost;
    long q_delmin = 0;
    if (b != NULL) {
        b->delcost += p->delcost;
        q_delmin = max(q_delmin, b->delmin - b->delcost);
    }
    if (c != NULL) {
        q_delmin = max(q_delmin, c->delmin - c->delcost);
    }
    long p_delmin = max(0, q_delmin - q_delcost);
    if (a != NULL) {
        p_delmin = max(p_delmin, a->delmin - a->delcost);
    }
    p->delcost = p_delcost;
    p->delmin = p_delmin;
    q->delcost = q_delcost;
    q->delmin = q_delmin;


    if ((q->l=b) != NULL) {
        b->p = q;
    }
    p->r = q;
    q->p = p;
    if ((p->p=r) != NULL) {
        if (r->l == q) r->l = p;
        else if (r->r == q) r->r = p;
    }
}


/*         q                   p
 *       /  \                /  \
 *     a     p     ->      q     c
 *         /  \          /  \
 *       b     c       a     b
 */
void rotL (dynamic_node_t* p) {
    dynamic_node_t* q = p->p;
    dynamic_node_t* r = q->p;
    dynamic_node_t* a = q->l;
    dynamic_node_t* b = p->l;
    dynamic_node_t* c = p->r;
    long p_delcost = p->delcost + q->delcost;
    long q_delcost = -p->delcost;
    long q_delmin = 0;
    if (b != NULL) {
        b->delcost += p->delcost;
        q_delmin = max(q_delmin, b->delmin - b->delcost);
    }
    if (a != NULL) {
        q_delmin = max(q_delmin, a->delmin - a->delcost);
    }
    long p_delmin = max(0, q_delmin - q_delcost);
    if (c != NULL) {
        p_delmin = max(p_delmin, c->delmin - c->delcost);
    }
    p->delcost = p_delcost;
    p->delmin = p_delmin;
    q->delcost = q_delcost;
    q->delmin = q_delmin;


    if ((q->r=b) != NULL) {
        q->r->p = q;
    }
    p->l = q;
    q->p = p;
    if ((p->p=r) != NULL) {
        if (r->l == q) r->l = p;
        else if (r->r == q) r->r = p;
    }
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
}

/* This makes node q the root of the virtual tree, and also q is the
   leftmost node in its splay tree*/
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

        dynamic_node_t* t = p->r;
        if (t != NULL) {
            p->delmin = max(p->delmin, t->delmin - t->delcost);
        }

        r = p;
    }

    splay(q);
}

/* assuming p and q are nodes in different trees and
   that p is a root of its tree, this links p to q */
void link(long pid, long qid, arc* edge) {
    dynamic_node_t* p = &dTree.nodes[pid];
    dynamic_node_t* q = &dTree.nodes[qid];

    expose(p);
    if (p->r != NULL) {
        // p is not a root. Error
        return;
    }

    dynamic_node_t* a = p->l;
    if (a != NULL) {
        a->delcost += p->delcost;
    }
    p->edge = edge;
    p->delcost = edge->resCap;
    p->delmin = 0;
    if (a != NULL) {
        a->delcost -= p->delcost;
        p->delmin = max(p->delmin, a->delmin - a->delcost);
    }

    p->p = q;

    // find the root
    dTree.cur_node = root(q->id);
}


/* this returns the id of the node that is the root of the tree containing p */
long root(long pid) {
    dynamic_node_t* p = &dTree.nodes[pid];
    expose(p);
    while(p->r != NULL) p = p->r;
    splay(p);
    return p->id;
}










// TODO: use splay inside before and after, not sure if it's needed
long before(long vid) {
    dynamic_node_t* v = &dTree.nodes[vid];
    if (v->l != NULL) {
        dynamic_node_t* u = v->l;
        for (;u->r != NULL;) {
            u = u->r;
        }
        return u->id;
    }
    else if (v->p != NULL && v->p->r == v) {
        return v->p->id;
    }
    else {
        return -1; // v is the head of the path, the leftmost node of the splay tree
    }
}



long after(long vid) {
    dynamic_node_t* v = &dTree.nodes[vid];
    if (v->r != NULL) {
        dynamic_node_t* u = v->r;
        for (;u->l != NULL;) {
            u = u->l;
        }
        return u->id;
    }
    else if (v->p != NULL && v->p->l == v) {
        return v->p->id;
    }
    else {
        return -1; // v is the tail of the path, the rightmost node of the splay tree
    }
}


// long pMinCost(dynamic_path_t* p) {
//     return nCost(p->root) - p->root->delmin;
// }


long nMinCost(long vid, dynamic_node_t** rootptr) {
    dynamic_node_t* v = &dTree.nodes[vid];
    long rootid = root(vid);
    dynamic_node_t* r = &dTree.nodes[rootid];
    if (rootptr != NULL) {
        *rootptr = r;
    }
    if (v == r) {
        // the path contains a single node
        return inf;
    }
    splay(r);
    // expose(v);
    dynamic_node_t* l = r->l;
    return r->delcost + l->delcost - l->delmin;
    // return v->delcost - v->delmin;
}


long nCost(long vid) {
    dynamic_node_t* v = &dTree.nodes[vid];
    splay(v);
    return v->delcost;
}


void pUpdate(long pid, long x) {
    dynamic_node_t* p = &dTree.nodes[pid];
    long rootid = root(pid);
    dynamic_node_t* r = &dTree.nodes[rootid];
    splay(r);
    dynamic_node_t* l = r->l;
    if (l != NULL) {
        l->delcost += x;
    }
}


void cut(long vid) {
    dynamic_node_t* v = &dTree.nodes[vid];
    expose(v);

    dynamic_node_t* r = v->r;
    dynamic_node_t* l = v->l;
    if (r != NULL) {
        r->p = NULL;
        v->r = NULL;
        r->delcost += v->delcost;
    }

    v->edge->resCap = v->delcost;

    if (l != NULL) {
        l->delcost += v->delcost;
    }
    v->delcost = v->delmin = 0;

    dTree.cur_node = vid;
}


void cutEdge(long vid) {
    dynamic_node_t* v = &dTree.nodes[vid];
    dynamic_node_t* u;
    while(nMinCost(vid, &u) == 0) {
        if (u == v) {
            //error
            break;
        }
        dynamic_node_t* w;
        long cost = u->delcost;
        u = u->l;
        cost += u->delcost;
        while(true) {
            if ((w = u->r) != NULL && cost + w->delcost - w->delmin == 0) {
                u = w;
                cost += u->delcost;
            }
            else if (cost == 0) {
                break;
            }
            else if ((w = u->l) != NULL && cost + w->delcost - w->delmin == 0) {
                u = w;
                cost += u->delcost;
            }
            // printf("while loop, uid = %ld, cost = %ld\n", u->id, cost);
        }
        // printf("cut = %ld\n", u->id);
        cut(u->id);
    }
    //dTree.cur_node = root(vid);
}


// store the path from p to the root
void findPath(long pid, long* a, long* b, long* cost) {
    dynamic_node_t* p = &dTree.nodes[pid];
    long pcost = nMinCost(pid, NULL);
    *b = before(root(pid));
    *a = after(pid);
    *cost = pcost;
    // printf("cost = %ld\n", pcost);
    pUpdate(pid, -pcost);
    // printf("update done\n");
    cutEdge(pid);
}
