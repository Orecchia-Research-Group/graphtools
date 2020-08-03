#include "types.h"


typedef struct {
//naive
    long id;
    long cost;// replace this
    arc* edge;
    node* node;
    void* prev_node;
    void* nex_node;
    // void* start_edge; // first adjacent edge in the adjacency list
    // void* cur_edge;
} dt_node_t;

typedef struct {
    dt_node_t* from;
    dt_node_t* to;
    long cost;
} dt_edge_t;

typedef struct {
//naive
    dt_node_t* head;
    dt_node_t* tail;
} dt_path_t;


typedef struct {
//naive
    long n,m;
    dt_node_t* nodes;
    dt_path_t** belongTo;

} dynamic_tree_t;


dynamic_tree_t dTree;







void init(long nodes, long edges);

void cleanUp();


long before(long v);
long after(long v);
dt_path_t* path(long v);
long head(dt_path_t* p);
long tail(dt_path_t* p);
long pCost(long v);
dt_path_t* concatenate(dt_path_t* p, dt_path_t* q, arc* edge);
long pMinCost(dt_path_t* p);
void pUpdate(dt_path_t* p, long x);
void split(dt_node_t* v);
void cut(dt_node_t* v);
void savePath(dt_node_t* a, dt_node_t* b, long cost);
void cutEdges(dt_path_t* p);
void findPath(dt_path_t* p, long* a, long* b, long* cost);

void splice(dt_path_t* p);
void expose(dt_node_t* v);


// dt_node_t* before(dt_node_t* v);
// dt_node_t* after(dt_node_t* v);
// dt_path_t* path(dt_node_t* v);
// dt_node_t* head(dt_path_t* p);
// dt_node_t* tail(dt_path_t* p);
// long pCost(dt_node_t* v);
// dt_path_t* concatenate(dt_path_t* p, dt_path_t* q, long cost);
// dt_node_t* pMinCost(dt_path_t* p);
// void pUpdate(dt_path_t* p, long x);
// void split(dt_node_t* v);
// void cut(dt_node_t* v);
// void savePath(dt_node_t* a, dt_node_t* b, long cost);
// void cut_edges(dt_path_t* p);
// void findPath()
//
// void splice(dt_path_t* p);
// void expose(dt_node_t* v);