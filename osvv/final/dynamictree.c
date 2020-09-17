#include "dynamictree.h"
#include "stdlib.h"
#include "stddef.h"
#include "stdio.h"

#define min(a,b) ((a) < (b)) ? (a) : (b)
#define max(a,b) ((a) > (b)) ? (a) : (b)

#define inf (0x1f1f1f1f1f1f1f1fL)


dynamic_tree_t* dec_init(long nodes, long start_node) {
    dynamic_tree_t* dTree = calloc(1, sizeof(dynamic_tree_t));
    dTree->sz = nodes;
    dTree->nodes = calloc(nodes, sizeof(dynamic_node_t));
    for (long i = 0; i < dTree->sz; i++) { // (TODO: int or long)
        initNode(&dTree->nodes[i], i);
    }
    dTree->cur_node = start_node;
    return dTree;
}

void initNode(dynamic_node_t* node, long i) {
    node->id = i;
    node->left = node->right = node->parent = NULL;
    node->delmin = 0;
    node->delcost = inf;
    node->edge = NULL;
}

void cleanUp(dynamic_tree_t* dTree) {
    free(dTree->nodes);
    free(dTree);
}

bool isroot(dynamic_node_t* node) {
    return (node->parent==NULL) || (node->parent->left != node && node->parent->right != node);
}





/*             q              p
 *           /  \           /  \
 *         p     c   ->   a     q
 *       /  \                 /  \
 *     a     b              b     c
 */
void rotR (dynamic_node_t* p) {
    dynamic_node_t* q = p->parent;
    dynamic_node_t* r = q->parent;
    dynamic_node_t* a = p->left;
    dynamic_node_t* b = p->right;
    dynamic_node_t* c = q->right;
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


    if ((q->left=b) != NULL) {
        b->parent = q;
    }
    p->right = q;
    q->parent = p;
    if ((p->parent=r) != NULL) {
        if (r->left == q) r->left = p;
        else if (r->right == q) r->right = p;
    }
}


/*         q                   p
 *       /  \                /  \
 *     a     p     ->      q     c
 *         /  \          /  \
 *       b     c       a     b
 */
void rotL (dynamic_node_t* p) {
    dynamic_node_t* q = p->parent;
    dynamic_node_t* r = q->parent;
    dynamic_node_t* a = q->left;
    dynamic_node_t* b = p->left;
    dynamic_node_t* c = p->right;
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


    if ((q->right=b) != NULL) {
        q->right->parent = q;
    }
    p->left = q;
    q->parent = p;
    if ((p->parent=r) != NULL) {
        if (r->left == q) r->left = p;
        else if (r->right == q) r->right = p;
    }
}

void splay(dynamic_node_t* p) {
    while (!isroot(p)) {
        dynamic_node_t* q = p->parent;
        if (isroot(q)) {
            if (q->left == p)
                rotR(p);
            else
                rotL(p);
        }
        else {
            dynamic_node_t* r = q->parent;
            if (r->left == q) {
                if (q->left == p) {
                    rotR(q);
                    rotR(p);
                }
                else {
                    rotL(p);
                    rotR(p);
                }
            }
            else {
                if (q->right == p) {
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
    for (dynamic_node_t* p=q; p != NULL; p=p->parent) {
        splay(p);

        dynamic_node_t* s = p->left;
        if (s != NULL) {
            s->delcost += p->delcost;
        }

        p->left = r;


        p->delmin = 0;
        if (r != NULL) {
            r->delcost -= p->delcost;
            p->delmin = max(p->delmin, r->delmin - r->delcost);
        }

        dynamic_node_t* t = p->right;
        if (t != NULL) {
            p->delmin = max(p->delmin, t->delmin - t->delcost);
        }

        r = p;
    }

    splay(q);
}

/* assuming p and q are nodes in different trees and
   that p is a root of its tree, this links p to q */
void link(dynamic_tree_t* dTree, long pid, long qid, arc* edge) {
    dynamic_node_t* p = &dTree->nodes[pid];
    dynamic_node_t* q = &dTree->nodes[qid];

    expose(p);
    if (p->right != NULL) {
        // p is not a root. Error
        return;
    }

    dynamic_node_t* a = p->left;
    if (a != NULL) {
        a->delcost += p->delcost;
    }
    p->edge = edge;
    p->delcost = edge->cap - edge->resCap;
    p->delmin = 0;
    if (a != NULL) {
        a->delcost -= p->delcost;
        p->delmin = max(p->delmin, a->delmin - a->delcost);
    }

    p->parent = q;

    // find the root
    dTree->cur_node = root(dTree, q->id);
}


/* this returns the id of the node that is the root of the tree containing p */
long root(dynamic_tree_t* dTree, long pid) {
    dynamic_node_t* p = &dTree->nodes[pid];
    expose(p);
    while(p->right != NULL) p = p->right;
    splay(p);
    return p->id;
}










// TODO: use splay inside before and after, not sure if it's needed
long before(dynamic_tree_t* dTree, long vid) {
    dynamic_node_t* v = &dTree->nodes[vid];
    if (v->left != NULL) {
        dynamic_node_t* u = v->left;
        for (;u->right != NULL;) {
            u = u->right;
        }
        return u->id;
    }
    else if (v->parent != NULL && v->parent->right == v) {
        return v->parent->id;
    }
    else {
        return -1; // v is the head of the path, the leftmost node of the splay tree
    }
}



long after(dynamic_tree_t* dTree, long vid) {
    dynamic_node_t* v = &dTree->nodes[vid];
    if (v->right != NULL) {
        dynamic_node_t* u = v->right;
        for (;u->left != NULL;) {
            u = u->left;
        }
        return u->id;
    }
    else if (v->parent != NULL && v->parent->left == v) {
        return v->parent->id;
    }
    else {
        return -1; // v is the tail of the path, the rightmost node of the splay tree
    }
}


// long pMinCost(dynamic_path_t* p) {
//     return nCost(p->root) - p->root->delmin;
// }


long nMinCost(dynamic_tree_t* dTree, long vid) {
    dynamic_node_t* v = &dTree->nodes[vid];
    expose(v);
    return v->delcost - v->delmin;
    // return v->delcost - v->delmin;
}


long nCost(dynamic_tree_t* dTree, long vid) {
    dynamic_node_t* v = &dTree->nodes[vid];
    splay(v);
    return v->delcost;
}


void pUpdate(dynamic_tree_t* dTree, long pid, long x) {
    dynamic_node_t* p = &dTree->nodes[pid];
    expose(p);
    p->delcost += x;
}


void cut(dynamic_tree_t* dTree, long vid) {
    dynamic_node_t* v = &dTree->nodes[vid];
    expose(v);

    dynamic_node_t* r = v->right;
    dynamic_node_t* l = v->left;
    if (r != NULL) {
        r->parent = NULL;
        v->right = NULL;
        r->delcost += v->delcost;
    }

    v->edge->resCap = v->edge->cap - v->delcost;

    v->delmin = 0;
    if (l != NULL) {
        l->delcost += v->delcost;
        l->delcost -= inf;
        v->delmin = max(v->delmin, l->delmin - l->delcost);
    }
    v->delcost = inf;

    dTree->cur_node = vid;
}


void cutEdge(dynamic_tree_t* dTree, long vid) {
    dynamic_node_t* v = &dTree->nodes[vid];
    while(nMinCost(dTree, vid) == 0) {
        dynamic_node_t* u = v;
        dynamic_node_t* w;
        long cost = u->delcost;
        while(true) {
            if ((w = u->right) != NULL && cost + w->delcost - w->delmin == 0) {
                u = w;
                cost += u->delcost;
            }
            else if (cost == 0) {
                break;
            }
            else if ((w = u->left) != NULL && cost + w->delcost - w->delmin == 0) {
                u = w;
                cost += u->delcost;
            }
        }
        cut(dTree, u->id);
    }
}


// store the path from p to the root
void findPath(dynamic_tree_t* dTree, long pid, long* a, long* b, long* cost) {
    long pcost = nMinCost(dTree, pid);
    *b = before(dTree, root(dTree, pid));
    *a = after(dTree, pid);
    *cost = pcost;
    pUpdate(dTree, pid, -pcost);
    cutEdge(dTree, pid);
}
