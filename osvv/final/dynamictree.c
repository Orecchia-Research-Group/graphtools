#include "dynamictree.h"
#include "stdlib.h"
#include "stddef.h"
#include "stdio.h"

#define min(a,b) ((a) < (b)) ? (a) : (b)
#define max(a,b) ((a) > (b)) ? (a) : (b)

#define inf (0x1f1f1f1f1f1f1f1fL)

void dt_print_op_stat(dynamic_tree_t* dTree) {
    long tot_rot_cnt = 0;
    long tot_splay_cnt = 0;
    long tot_expose_cnt = 0;
    for (int i = 0; i < dTree->sz; i++){
        tot_rot_cnt += dTree->d_nodes[i].rot_cnt;
        tot_splay_cnt += dTree->d_nodes[i].splay_cnt;
        tot_expose_cnt += dTree->d_nodes[i].expose_cnt;
    }
    fprintf(stderr, "\n\n\n------------------------------------\n");
    fprintf(stderr, "Size of the dynamic tree = %ld\n", dTree->sz);
    fprintf(stderr, "Number of links = %ld\n", dTree->link_cnt);
    fprintf(stderr, "Number of cuts = %ld\n", dTree->cut_cnt);
    fprintf(stderr, "Number of rotations = %ld\n", tot_rot_cnt);
    fprintf(stderr, "Number of splay = %ld\n", tot_splay_cnt);
    fprintf(stderr, "Number of expose = %ld\n", tot_expose_cnt);
    fprintf(stderr, "\n------------------------------------\n\n\n");
}

dynamic_tree_t* dt_init(long num_nodes, node* nodes, node* start_node) {
    dynamic_tree_t* dTree = calloc(1, sizeof(dynamic_tree_t));
    dTree->sz = num_nodes;
    dTree->d_nodes = calloc(num_nodes, sizeof(dynamic_node_t));
    for (long i = 0; i < dTree->sz; i++) { // (TODO: int or long)
        dt_initNode(&dTree->d_nodes[i], i);
    }
    dTree->nodes = nodes;
    dTree->source = start_node;
    dTree->d_source = dt_to_d_node(dTree, dTree->source);
    dTree->cur_node = start_node;
    dTree->d_cur_node = dt_to_d_node(dTree, dTree->cur_node);

    dTree->link_cnt = 0;
    dTree->cut_cnt = 0;
    return dTree;
}

void dt_initNode(dynamic_node_t* node, long i) {
    node->id = i;
    node->left = node->right = node->parent = NULL;
    node->delmin = 0;
    node->delcost = inf;
    node->edge = NULL;

    node->rot_cnt = 0;
    node->splay_cnt = 0;
    node->expose_cnt = 0;
}

void dt_cleanUp_d_node(dynamic_node_t* d_node, bool* flag, long parent_cost) {
    flag[d_node->id] = true;
    long cost = parent_cost + d_node->delcost;
    if (d_node->edge != NULL) {
        d_node->edge->resCap = d_node->edge->cap - cost;
    }
    if (d_node->left != NULL) {
        dt_cleanUp_d_node(d_node->left, flag, cost);
    }
    if (d_node->right != NULL) {
        dt_cleanUp_d_node(d_node->right, flag, cost);
    }
}

void dt_cleanUp(dynamic_tree_t* dTree) {
    bool* flag = calloc(dTree->sz, sizeof(bool));
    for (long i = 0; i < dTree->sz; i++) {
        if (flag[i] == false) {
            dt_splay(&dTree->d_nodes[i]);
            dt_cleanUp_d_node(&dTree->d_nodes[i], flag, 0);
        }
    }
    free(flag);

    free(dTree->d_nodes);
    free(dTree);
}

dynamic_node_t* dt_to_d_node(dynamic_tree_t* dTree, node* p) {
    // fprintf(stderr, "nodes = %p, d_nodes = %p, node p = %p\n", dTree->nodes, dTree->d_nodes, p);
    if (p == NULL) return NULL;
    return dTree->d_nodes + (p - dTree->nodes);
}

node* dt_to_node(dynamic_tree_t* dTree, dynamic_node_t* p) {
    // fprintf(stderr, "nodes = %p, d_nodes = %p, d_node p = %p\n", dTree->nodes, dTree->d_nodes, p);
    if (p == NULL) return NULL;
    return dTree->nodes + (p - dTree->d_nodes);
}

bool dt_isroot(dynamic_node_t* node) {
    return (node->parent==NULL) || (node->parent->left != node && node->parent->right != node);
}

/*             q              p
 *           /  \           /  \
 *         p     c   ->   a     q
 *       /  \                 /  \
 *     a     b              b     c
 */
void dt_rotR (dynamic_node_t* p) {
    p->rot_cnt++;

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
void dt_rotL (dynamic_node_t* p) {
    p->rot_cnt++;

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

void dt_splay(dynamic_node_t* p) {
    p->splay_cnt++;
    while (!dt_isroot(p)) {
        dynamic_node_t* q = p->parent;
        if (dt_isroot(q)) {
            if (q->left == p)
                dt_rotR(p);
            else
                dt_rotL(p);
        }
        else {
            dynamic_node_t* r = q->parent;
            if (r->left == q) {
                if (q->left == p) {
                    dt_rotR(q);
                    dt_rotR(p);
                }
                else {
                    dt_rotL(p);
                    dt_rotR(p);
                }
            }
            else {
                if (q->right == p) {
                    dt_rotL(q);
                    dt_rotL(p);
                }
                else {
                    dt_rotR(p);
                    dt_rotL(p);
                }
            }
        }
    }
}

/* This makes node q the root of the virtual tree, and also q is the
   leftmost node in its splay tree*/
void dt_expose(dynamic_node_t* q) {
    q->expose_cnt++;
    dynamic_node_t* r = NULL;
    for (dynamic_node_t* p=q; p != NULL; p=p->parent) {
        dt_splay(p);

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

    dt_splay(q);
}

/* assuming p and q are nodes in different trees and
   that p is a root of its tree, this links p to q */
void dt_d_link(dynamic_tree_t* dTree, dynamic_node_t* p, dynamic_node_t* q, arc* edge) {
    dTree->link_cnt++;
    dt_expose(p);
    if (p->right != NULL) {
        // p is not a root. Error
        long pid = p - dTree->d_nodes;
        long qid = q - dTree->d_nodes;
        fprintf(stderr, "Linking from node %ld to node %ld failed, because node %ld is not a root\n", pid, qid, pid);
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
    dTree->d_cur_node = dt_d_root(q);
    dTree->cur_node = dt_to_node(dTree, dTree->d_cur_node);
}

void dt_link(dynamic_tree_t* dTree, node* p, node* q, arc* edge) {
    dt_d_link(dTree, dt_to_d_node(dTree, p), dt_to_d_node(dTree, q), edge);
}


/* this returns the id of the node that is the root of the tree containing p */
dynamic_node_t* dt_d_root(dynamic_node_t* p) {
    dt_expose(p);
    while(p->right != NULL) p = p->right;
    dt_splay(p);
    return p;
}

node* dt_root(dynamic_tree_t* dTree, node* p) {
    dynamic_node_t* d_p = dt_to_d_node(dTree, p);
    d_p = dt_d_root(d_p);
    p = dt_to_node(dTree, d_p);
    return p;
}

// TODO: use dt_splay inside before and after, not sure if it's needed
dynamic_node_t* dt_d_before(dynamic_node_t* p) {
    if (p->left != NULL) {
        dynamic_node_t* u = p->left;
        for (;u->right != NULL;) {
            u = u->right;
        }
        return u;
    }
    else if (p->parent != NULL && p->parent->right == p) {
        return p->parent;
    }
    else {
        // fprintf(stderr, "node %p is the head, before(head) == NULL\n", p);
        return NULL; // p is the head of the path, the leftmost node of the dt_splay tree
    }
}

dynamic_node_t* dt_d_after(dynamic_node_t* p) {
    if (p->right != NULL) {
        dynamic_node_t* u = p->right;
        for (;u->left != NULL;) {
            u = u->left;
        }
        return u;
    }
    else if (p->parent != NULL && p->parent->left == p) {
        return p->parent;
    }
    else {
        // fprintf(stderr, "node %p is the tail, after(tail) == NULL\n", p);
        return NULL; // p is the tail of the path, the rightmost node of the dt_splay tree
    }
}

node* dt_before(dynamic_tree_t* dTree, node* p) {
    dt_expose(dTree->d_source);
    dynamic_node_t* d_p = dt_to_d_node(dTree, p);
    d_p = dt_d_before(d_p);
    p = dt_to_node(dTree, d_p);
    return p;
}

node* dt_after(dynamic_tree_t* dTree, node* p) {
    dt_expose(dTree->d_source);
    dynamic_node_t* d_p = dt_to_d_node(dTree, p);
    d_p = dt_d_after(d_p);
    p = dt_to_node(dTree, d_p);
    return p;
}

long dt_nMinCost(dynamic_node_t* p) {
    dt_expose(p);
    return p->delcost - p->delmin;
}

long dt_nCost(dynamic_node_t* p) {
    dt_splay(p);
    return p->delcost;
}

void dt_pUpdate(dynamic_node_t* p, long x) {
    dt_expose(p);
    p->delcost += x;
}

void dt_d_cut(dynamic_tree_t* dTree, dynamic_node_t* p) {
    dTree->cut_cnt++;
    dt_expose(p);

    dynamic_node_t* r = p->right;
    dynamic_node_t* l = p->left;
    if (r != NULL) {
        r->parent = NULL;
        p->right = NULL;
        r->delcost += p->delcost;
    }

    p->edge->resCap = p->edge->cap - p->delcost;
    p->edge = NULL;

    p->delmin = 0;
    if (l != NULL) {
        l->delcost += p->delcost;
        l->delcost -= inf;
        p->delmin = max(p->delmin, l->delmin - l->delcost);
    }
    p->delcost = inf;

    dTree->d_cur_node = p;
    dTree->cur_node = dt_to_node(dTree, dTree->d_cur_node);
}

void dt_cut(dynamic_tree_t* dTree, node* p) {
    dynamic_node_t* d_p = dt_to_d_node(dTree, p);
    dt_d_cut(dTree, d_p);
}

void dt_d_cutEdge(dynamic_tree_t* dTree, dynamic_node_t* p) {
    while(dt_nMinCost(p) == 0) {
        dynamic_node_t* u = p;
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
        dt_d_cut(dTree, u);
    }
}

void dt_cutEdge(dynamic_tree_t* dTree, node* p) {
    dynamic_node_t* d_p = dt_to_d_node(dTree, p);
    dt_d_cutEdge(dTree, d_p);
}

// store the path from p to the root
void dt_findPath(dynamic_tree_t* dTree, node** a, node** b, long* cost) {
    node* p = dTree->source;
    long pcost = dt_nMinCost(dt_to_d_node(dTree, p));
    *b = dt_before(dTree, dt_root(dTree, p));
    *a = dt_after(dTree, p);
    *cost = pcost;
    dt_pUpdate(dt_to_d_node(dTree, p), -pcost);
    dt_cutEdge(dTree, p);
}
