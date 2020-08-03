#include "dectree.h"
#include "stdlib.h"
#include "stddef.h"
#include "stdbool.h"


#define min(a,b) ((a) < (b))? (a) : (b)


void init(long nodes, long edges) {
    dTree.n = nodes;
    dTree.m = edges;
    dTree.nodes = (dt_node_t*) calloc(dTree.n + 2, sizeof(dt_node_t));
    dTree.belongTo = (dt_path_t**) calloc(dTree.n + 2, sizeof(dt_path_t*));

    for (int i = 0; i < dTree.n + 2; i++) {
        // dTree.nodes[i] = calloc(1, sizeof(dt_node_t));
        dTree.nodes[i].id = i;
        dTree.nodes[i].edge = NULL;
        dTree.nodes[i].node = NULL; //TODO: link to nodes in flow.c
        dTree.nodes[i].prev_node = NULL;
        dTree.nodes[i].next_node = NULL;
        //TODO: implement adjacency list through start_edge, cur_edge


        dTree.belongTo[i] = calloc(1, sizeof(dt_path_t));
        dTree.belongTo[i]->head = dTree.belongTo[i]-> tail = &dTree.nodes[i];
    }
}


void cleanUp() {

}


// TODO: add error and warning messages
long before(long v) {
    return dTree.nodes[v].prev_node->id;
}

long after(long v) {
    return dTree.nodes[v].next_node->id;
}

dt_path_t* path(long v) {
    return dTree.belongTo[v];
}

long head(dt_path_t* p) {
    return (p != NULL)? p->head->id : -1;
}

long tail(dt_path_t* p) {
    return (p != NULL)? p->tail->id : -1;
}

long pCost(long v) {
    if(dTree.nodes[v].edge == NULL) return -1;
    return dTree.nodes[v].edge->resCap;
}

dt_path_t* concatenate(dt_path_t* p, dt_path_t* q, arc* edge) {
    dt_node_t* p_tail = p->tail;
    dt_node_t* q_head = q->head;
    dt_node_t* q_tail = q->tail;
    p_tail->next_node = q_head;
    p_tail->edge = edge;
    q_head->prev_node = p_tail;
    dt_node_t* q_node = q_head;
    while (q_node != NULL) {
        dTree.belongTo[q_node->id] = p;
        q_node = q_node->next_node;
    }
    p->tail = q_tail;
    free(q);
    //TODO, add cost information somewhere
    return NULL;
}

long pMinCost(dt_path_t* p) {
    long min_cost = 0x3f3f3f3f3f3f3f3f;
    dt_node_t* min_node = NULL;
    dt_node_t* node = p->head;
    dt_node_t* tail = p->tail;
    while(node != tail) {
        if (min_cost >= node->edge->resCap) {
            min_cost = node->edge->resCap;
            min_node = node;
        }
        node = node->next_node;
    }
    return min_node->id;
}

void pUpdate(dt_path_t* p, long x) {
    dt_node_t* node = p->head;
    dt_node_t* tail = p->tail;
    while(node != tail) {
        node->edge->resCap += x;
        node = node->next_node;
    }
}

void split(dt_node_t* v) {
    //not used
}

void cut(dt_node_t* v) {
    dt_node_t* w = v->next_node;
    dt_path_t* v_path = dTree.belongTo[v->id];
    dt_path_t* w_path = calloc(1, sizeof(dt_path_t));
    dt_node_t* w_tail = v_path->tail;
    v->next_node = NULL;
    w->prev_node = NULL;
    v_path->tail = v;

    w_path->head = w;
    w_path->tail = w_tail;
    dt_node_t* node = w;
    while(node != NULL) {
        dTree.belongTo[node->id] = w_path;
        node = node->next_node;
    }
}

void savePath(dt_node_t* a, dt_node_t* b, long cost) {

}

void cutEdges(dt_path_t* p) {
    while(true) {
        long node = pMinCost(p);
        if (node == p->tail->id) {
            break;
        }
        if (node && dTree.nodes[node]->edge->resCap == 0) {
            cut(&dTree.nodes[node]);
        }
        else {
            break;
        }
    }
}


void findPath(dt_path_t* p, long* a, long* b, long* cost) {
    *a = after(p->head->id);
    *b = before(p->tail->id);
    long min_node = pMinCost(p);
    long flow = dTree.nodes[min_node].edge->resCap;
    *cost = flow;
    pUpdate(p, -flow);
    cutEdges(p);
}



void splice(dt_path_t* p) {
    //not used
}

void expose(dt_node_t* v) {
    //not used
}
