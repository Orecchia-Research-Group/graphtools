
typedef struct {
    long id;
} path_t;


typedef struct {
    
} dynamic_tree_t;



long before(long v);
long after(long v);
path_t* path(long v);
long head(path_t* p);
long tail(path_t* p);
long pCost(long v);
path_t* concatenate(path_t* p, path_t* q, long cost);
long pMinCost(path_t* p);
void pUpdate(path_t* p, long x);
void split(long v);
void savePath(long a, long b, long cost);

void splice(path_t* p);
void expose(long v);
