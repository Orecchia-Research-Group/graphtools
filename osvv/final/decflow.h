
typedef struct {
//naive
    long id;
    node* head;
} path_t;


typedef struct {
//naive
    long n,m;
    path_t** belongTo;

} dynamic_tree_t;


dynamic_tree_t dTree;







void init(long nodes, long edges);

void cleanUp();



node* before(node* v);
node* after(node* v);
path_t* path(node* v);
node* head(path_t* p);
node* tail(path_t* p);
long pCost(node* v);
path_t* concatenate(path_t* p, path_t* q, long cost);
long pMinCost(path_t* p);
void pUpdate(path_t* p, long x);
void split(node* v);
void savePath(node* a, node* b, long cost);

void splice(path_t* p);
void expose(node* v);
