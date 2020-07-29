#include "decflow.h"


#define min(a,b) ((a) < (b))? (a) : (b)


void init(long nodes, long edges) {
    dTree.n = nodes;
    dTree.m = edges;
    dTree.nodes = (dt_node_t**) calloc(dTree.n + 2, sizeof(dt_node_t*));
    dTree.belongTo = (dt_path_t**) calloc(dTree.n + 2, sizeof(dt_path_t*));

    for (int i = 0; i < dTree.n + 2; i++) {
        dTree.nodes[i] = calloc(1, sizeof(dt_node_t*));
        dTree.nodes[i]->id = i;
        dTree.nodes[i]->cost = -1;
        dTree.nodes[i]->node = NULL; //TODO: link to nodes in flow.c
        dTree.nodes[i]->prev_node = NULL;
        dTree.nodes[i]->nex_node = NULL;
        //TODO: implement adjacency list through start_edge, cur_edge


        dTree.belongTo[i] = calloc(1, sizeof(dt_path_t*));
        dTree.belongTo[i]->head = dTree.belongTo[i]-> tail = dTree.nodes[i];
    }
}


void cleanUp() {

}


// TODO: add error and warning messages
dt_node_t* before(dt_node_t* v) {
    return (v != NULL)? v->prev_node : NULL;
}

dt_node_t* after(dt_node_t* v) {
    return (v != NULL)? v->nex_node : NULL;
}

dt_path_t* path(dt_node_t* v) {
    return (v != NULL)? dTree.belongTo[v->id] : NULL;
}

dt_node_t* head(dt_path_t* p) {
    return (p != NULL)? p->head : NULL;
}

dt_node_t* tail(dt_path_t* p) {
    return (p != NULL)? p->tail : NULL;
}

long pCost(dt_node_t* v) {
    return (v != NULL)? v->cost : -1;
}

dt_path_t* concatenate(dt_path_t* p, dt_path_t* q, long cost) {
    dt_node_t* p_tail = p->tail;
    dt_node_t* q_head = q->head;
    dt_node_t* q_tail = q->tail;
    p_tail->nex_node = q_head;
    p_tail->cost = cost;
    q_head->prev_node = p_tail;
    dt_node_t* q_node = q_head;
    while (q_node != NULL) {
        dTree.belongTo[q_node->id] = p;
        q_node = q_node->nex_node;
    }
    p->tail = q_tail;
    free(q);
    //TODO, add cost information somewhere
    return NULL;
}

dt_node_t* pMinCost(dt_path_t* p) {
    long min_cost = 0x3f3f3f3f3f3f3f3f;
    dt_node_t* min_node = NULL;
    dt_node_t* node = p->head;
    dt_node_t* tail = p->tail;
    while(node != tail) {
        if (min_cost >= node->cost) {
            min_cost = node->cost;
            min_node = node;
        }
        node = node->nex_node;
    }
    return min_node;
}

void pUpdate(dt_path_t* p, long x) {
    dt_node_t* node = p->head;
    dt_node_t* tail = p->tail;
    while(node != tail) {
        node->cost += x;
        node = node->nex_node;
    }
}

void split(dt_node_t* v) {
    //not used
}

void cut(dt_node_t* v) {
    dt_node_t* w = v->nex_node;
    dt_path_t* v_path = dTree.belongTo[v->id];
    dt_path_t* w_path = calloc(1, sizeof(dt_path_t));
    dt_node_t* w_tail = v_path->tail;
    v->nex_node = NULL;
    w->prev_node = NULL;
    v_path->tail = v;

    w_path->head = w;
    w_path->tail = w_tail;
    dt_node_t* node = w;
    while(node != NULL) {
        dTree.belongTo[node->id] = w_path;
        node = node->nex_node;
    }
}

void savePath(dt_node_t* a, dt_node_t* b, long cost) {

}

void cut_edges(dt_path_t* p) {
    while(true) {
        dt_node_t* node = pMinCost(p);
        if (node && node->cost == 0) {
            cut(node);
            last_node = node;
        }
        else {
            break;
        }
    }
}


void splice(dt_path_t* p) {
    //not used
}

void expose(dt_node_t* v) {
    //not used
}
