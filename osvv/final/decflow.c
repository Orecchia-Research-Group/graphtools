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
    return (v != NULL)? dTree.belongTo[v.id] : NULL;
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
    p_tail->nex_node = q_head;
    q_head->prev_node = p_tail;
    dt_node_t* q_node = q_head;
    while (q_node != NULL) {
        dTree.belongTo[q_node->id] = p;
        q_node = q_node->nex_node;
    }
    //TODO, add cost information somewhere
    return NULL;
}

long pMinCost(dt_path_t* p) {
    long min_cost = 0x3f3f3f3f3f3f3f3f;
    dt_node_t* node = p->head;
    dt_node_t* tail = p->tail;
    while(node != tail) {
        min_cost = min(min_cost, node->cost);
    }
    return min_cost;
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

}

void savePath(dt_node_t* a, dt_node_t* b, long cost) {

}


void splice(dt_path_t* p) {

}

void expose(dt_node_t* v) {

}
