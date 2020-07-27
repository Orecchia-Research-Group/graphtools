
typedef struct {
//naive
    int id;
    long cost;
    node* node;
    void* prev_node;
    void* nex_node;
    void* start_edge; // first adjacent edge in the adjacency list
    void* cur_edge;
} dt_node_t;

typedef struct {
    dt_node_t* from;
    dt_node_t* to;
    long cost;
    void* nex_edge;
} dt_edge_t;

typedef struct {
//naive
    dt_node_t* head;
    dt_node_t* tail;
} dt_path_t;


typedef struct {
//naive
    long n,m;
    dt_node_t** nodes;
    dt_path_t** belongTo;

} dynamic_tree_t;


dynamic_tree_t dTree;







void init(long nodes, long edges);

void cleanUp();



dt_node_t* before(dt_node_t* v);
dt_node_t* after(dt_node_t* v);
dt_path_t* path(dt_node_t* v);
dt_node_t* head(dt_path_t* p);
dt_node_t* tail(dt_path_t* p);
long pCost(dt_node_t* v);
dt_path_t* concatenate(dt_path_t* p, dt_path_t* q, long cost);
long pMinCost(dt_path_t* p);
void pUpdate(dt_path_t* p, long x);
void split(dt_node_t* v);
void savePath(dt_node_t* a, dt_node_t* b, long cost);

void splice(dt_path_t* p);
void expose(dt_node_t* v);
