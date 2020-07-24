#include "decflow.h"



void init(long nodes, long edges) {
    dTree.n = nodes;
    dTree.m = edges;
    dTree.belongTo = (path_t *) calloc(dTree.n + 2, sizeof(path_t *));

    for (int i = 0; i < dTree.n + 2; i++) {
        dTree.belongTo[i] = calloc(1, sizeof(path_t));
        dTree.belongTo[i]->id = i;
    }
}

void cleanUp() {

}






node* before(node* v) {
    return NULL;
}

node* after(node* v) {
    return NULL;
}

path_t* path(node* v) {
    return NULL;
}

node* head(path_t* p) {
    return NULL;
}

node* tail(path_t* p) {
    return NULL;
}

long pCost(node* v) {
    return -1;
}

path_t* concatenate(path_t* p, path_t* q, long cost) {
    return NULL;
}

long pMinCost(path_t* p) {
    return -1;
}

void pUpdate(path_t* p, long x) {

}

void split(node* v) {

}

void savePath(node* a, node* b, long cost) {

}


void splice(path_t* p) {

}

void expose(node* v) {

}
