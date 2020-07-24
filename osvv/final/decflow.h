
typedef struct {
    long id;
} path_t;


typedef struct {

} dynamic_tree_t;



long before(node* v);
long after(node* v);
path_t* path(node* v);
long head(path_t* p);
long tail(path_t* p);
long pCost(node* v);
path_t* concatenate(path_t* p, path_t* q, long cost);
long pMinCost(path_t* p);
void pUpdate(path_t* p, long x);
void split(node* v);
void savePath(node* a, node* b, long cost);

void splice(path_t* p);
void expose(node* v);
