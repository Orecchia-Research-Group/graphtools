/*
C function: hipr

PURPOSE: run maxflow computation on directed graph.

This function includes modifications by Reid Andersen and Satish Rao to output
the mincut and the demand graph routed in case of matching computations.

*/

/*
C function: hipr.

INPUTS: Note that vertex indices go from 1 to n.
 - n: number of vertices
 - m: number of arcs
 - tails: pointer to array of tails of the arcs. tails[i] is the tail of the (i-1)th arc.
 - heads: pointer to array of heads of the arcs.
 - weights: pointer to array of weights of the arcs.
 - s: index of source. (should be n-1)
 - t: index of sink . (should be n)
 - output_set: pointer to array of size n-2 (no sink and source) filled in by hipr to be mask for mincut.
 - mheads: pointer to array of heads of arcs of routed matching.
 - mtails: pointer to array of tails of arcs of routed matching.
 - mweights: pointer to array of weights of arcs of routed matching.
 - nedges: pointer to number of edges in mgraph
 - fflow: pointer to long which becomes equal to flow routed.
 - route_flag: set to 1 if we want to receive mgraph as output
*/




/* Maximum flow - highest lavel push-relabel algorithm */
/* COPYRIGHT C 1995, 2000 by IG Systems, Inc., igsys@eclipse.net */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
/*#include <values.h>*/
#define MAXLONG 1000000000

#include "types.h"          /* type definitions */
/* #include "flow.h" */
#include "timer.h"        /* timing routine */
#include "dynamictree.h"
/* MATLAB INTEGRATION LIBS */



#define GLOB_UPDT_FREQ 0.5
#define ALPHA 6
#define BETA 12

#define WHITE 0
#define GREY 1
#define BLACK 2

/* #define PRINT_CUT */
#define CHECK_SOLUTION

int loadflowproblem
        (
                long n,
                long m,
                long *tails,
                long *heads,
                long *weights,
                long s,
                long t,
                long *n_ad,
                long *m_ad,
                node **nodes_ad,
                arc **arcs_ad,
                long **cap_ad,
                node **source_ad,
                node **sink_ad,
                long *node_min_ad
        );


/* global variables */

long n;                    /* number of nodes */
long m;                    /* number of arcs */
long nm;                   /* n + ALPHA * m */
long nMin;                 /* smallest node id */
node *nodes;               /* array of nodes */
arc *arcs;                /* array of arcs */
bucket *buckets;             /* array of buckets */
cType *cap;                 /* array of capacities */
node *source;              /* source node pointer */
node *sink;                /* sink node pointer */
node   **queue;              /* queue for BFS */
node   **qHead, **qTail, **qLast;     /* queue pointers */
long dMax;                 /* maximum label */
long aMax;                 /* maximum active node label */
long aMin;                 /* minimum active node label */
double flow;                 /* flow value */
long pushCnt = 0;           /* number of pushes */
long relabelCnt = 0;       /* number of relabels */
long updateCnt = 0;       /* number of updates */
long gapCnt = 0;           /* number of gaps */
long gNodeCnt = 0;           /* number of nodes after gap */
float t1, t2;                 /* for saving times */
node *sentinelNode;        /* end of the node list marker */
arc *stopA;                  /* used in forAllArcs */
long workSinceUpdate = 0;      /* the number of arc scans since last update */
float globUpdtFreq;          /* global update frequency */

/* macros */
#define addedge(t, h, c)\
{\
           arc_first[t + 1] ++;\
           arc_first[h + 1] ++;\
           arc_tail[pos_current]        = t;\
           arc_tail[pos_current+1]      = h;\
           arc_current       -> head    = nodes + h;\
           arc_current       -> resCap    = c;\
           arc_current       -> rev  = arc_current + 1;\
           ( arc_current + 1 ) -> head    = nodes + t;\
           ( arc_current + 1 ) -> resCap  = 0;\
           ( arc_current + 1 ) -> rev  = arc_current;\
           arc_current += 2;\
           pos_current += 2;\
           j++;\
}

#define forAllNodes(i) for ( i = nodes; i != sentinelNode; i++ )
#define forAllArcs(i, a) for (a = (i)->first, stopA = (i+1)->first; a != stopA; a++)

#define nNode(i) ( (i) - nodes + nMin )
#define nArc(a)  ( ( a == NULL )? -1 : (a) - arcs )

#define min(a, b) ( ( (a) < (b) ) ? a : b )

/* FIFO queue for BFS macros */

#define qInit(n) \
{\
  qHead = qTail = queue;\
  qLast = queue + (n)-1;\
}

#define qEmpty ( qHead == qTail )

#define qEnqueue(i) \
{\
  *qTail = i;\
  if ( qTail == qLast ) qTail = queue;\
  else qTail++;\
}

#define qDequeue(i) \
{\
  i = qHead;\
  if ( qHead == qLast ) qHead = queue;\
  else qHead++;\
}


/*
   bucket macros:
   bucket's active node list is singly-linked
     operations aAdd, aRemove (from the front)
   bucket's inactive list is doubly-linked
     operations iAdd, iDelete (from arbitrary position)
*/

long i_dist;

#define aAdd(l, i)\
{\
  i->bNext = l->firstActive;\
  l->firstActive = i;\
  i_dist = i->d;\
  if (i_dist < aMin)\
    aMin = i_dist;\
  if (i_dist > aMax)\
    aMax = i_dist;\
  if (dMax < aMax)\
    dMax = aMax;\
}

/* i must be the first element */
#define aRemove(l, i)\
{\
  l->firstActive = i->bNext;\
}

node *i_next, *i_prev;
#define iAdd(l, i)\
{\
  i_next = l->firstInactive;\
  i->bNext = i_next;\
  i->bPrev = sentinelNode;\
  i_next->bPrev = i;\
  l->firstInactive = i;\
}

#define iDelete(l, i)\
{\
  i_next = i->bNext;\
  if (l->firstInactive == i) {\
    l->firstInactive = i_next;\
    i_next->bPrev = sentinelNode;\
  }\
  else {\
    i_prev = i->bPrev;\
    i_prev->bNext = i_next;\
    i_next->bPrev = i_prev;\
  }\
}

/* allocate datastructures, initialize related variables */

int allocDS() {

    nm = ALPHA * n + m;
    /*
    queue = (node**) calloc ( n, sizeof (node*) );
    if ( queue == NULL ) return ( 1 );
    qLast = queue + n - 1;
    qInit();
    */
    buckets = (bucket *) calloc(n + 2, sizeof(bucket));
    if (buckets == NULL) return (1);

    sentinelNode = nodes + n;
    sentinelNode->first = arcs + 2 * m;

    return (0);

} /* end of allocate */


void init() {
    node *i;        /* current node */
    int overflowDetected;
    bucket *l;
    arc *a;
#ifdef EXCESS_TYPE_LONG
    double testExcess;
#endif
#ifndef OLD_INIT
    unsigned long delta;
#endif

    /* initialize excesses */

    forAllNodes(i) {
        i->excess = 0;
        i->current = i->first;
        forAllArcs(i, a)a->resCap = cap[a - arcs];
    }

    for (l = buckets; l <= buckets + n - 1; l++) {
        l->firstActive = sentinelNode;
        l->firstInactive = sentinelNode;
    }

    overflowDetected = 0;
#ifdef EXCESS_TYPE_LONG
    testExcess = 0;
    forAllArcs(source,a) {
      if (a->head != source) {
        testExcess += a->resCap;
      }
    }
    if (testExcess > MAXLONG) {
      printf("c WARNING: excess overflow. See README for details.\nc\n");
      overflowDetected = 1;
    }
#endif
#ifdef OLD_INIT
    source -> excess = MAXLONG;
#else
    if (overflowDetected) {
        source->excess = MAXLONG;
    } else {
        source->excess = 0;
        forAllArcs(source, a) {
            if (a->head != source) {
                pushCnt++;
                delta = a->resCap;
                a->resCap -= delta;
                (a->rev)->resCap += delta;
                a->head->excess += delta;
            }
        }
    }

    /*  setup labels and buckets */
    l = buckets + 1;

    aMax = 0;
    aMin = n;

    forAllNodes(i) {
        if (i == sink) {
            i->d = 0;
            iAdd(buckets, i);
            continue;
        }
        if ((i == source) && (!overflowDetected)) {
            i->d = n;
        } else
            i->d = 1;
        if (i->excess > 0) {
            /* put into active list */
            aAdd(l, i);
        } else { /* i -> excess == 0 */
            /* put into inactive list */
            if (i->d < n) iAdd(l, i);
        }
    }
    dMax = 1;
#endif

    /*  dMax = n-1;
        flow = 0.0; */

} /* end of init */

void checkMax() {
    bucket *l;

    for (l = buckets + dMax + 1; l < buckets + n; l++) {
        assert(l->firstActive == sentinelNode);
        assert(l->firstInactive == sentinelNode);
    }
}

/* global update via backward breadth first search from the sink */

void globalUpdate() {

    node *i, *j;       /* node pointers */
    arc *a;           /* current arc pointers  */
    bucket *l, *jL;          /* bucket */
    long curDist, jD;
    long state;


    updateCnt++;

    /* initialization */

    forAllNodes(i)i->d = n;
    sink->d = 0;

    for (l = buckets; l <= buckets + dMax; l++) {
        l->firstActive = sentinelNode;
        l->firstInactive = sentinelNode;
    }

    dMax = aMax = 0;
    aMin = n;

    /* breadth first search */

    /* add sink to bucket zero */

    iAdd(buckets, sink);
    for (curDist = 0; 1; curDist++) {

        state = 0;
        l = buckets + curDist;
        jD = curDist + 1;
        jL = l + 1;
        /*
        jL -> firstActive   = sentinelNode;
        jL -> firstInactive  = sentinelNode;
        */

        if ((l->firstActive == sentinelNode) &&
            (l->firstInactive == sentinelNode))
            break;

        while (1) {

            switch (state) {
                case 0:
                    i = l->firstInactive;
                    state = 1;
                    break;
                case 1:
                    i = i->bNext;
                    break;
                case 2:
                    i = l->firstActive;
                    state = 3;
                    break;
                case 3:
                    i = i->bNext;
                    break;
                default:
                    assert(0);
                    break;
            }

            if (i == sentinelNode) {
                if (state == 1) {
                    state = 2;
                    continue;
                } else {
                    assert(state == 3);
                    break;
                }
            }

            /* scanning arcs incident to node i */
            forAllArcs(i, a) {
                if (a->rev->resCap > 0) {
                    j = a->head;
                    if (j->d == n) {
                        j->d = jD;
                        j->current = j->first;
                        if (jD > dMax) dMax = jD;

                        if (j->excess > 0) {
                            /* put into active list */
                            aAdd(jL, j);
                        } else {
                            /* put into inactive list */
                            iAdd(jL, j);
                        }
                    }
                }
            } /* node i is scanned */
        }
    }

} /* end of global update */

/* second stage -- preflow to flow */
void stageTwo()
/*
   do dsf in the reverse flow graph from nodes with excess
   cancel cycles if found
   return excess flow in topological order
*/

/*
   i->d is used for dfs labels
   i->bNext is used for topological order list
   buckets[i-nodes]->firstActive is used for DSF tree
*/

{
    node *i, *j, *tos, *bos, *restart, *r;
    arc *a;
    cType delta;

    /* deal with self-loops */
    forAllNodes(i) {
        forAllArcs(i, a)if (a->head == i) {
                a->resCap = cap[a - arcs];
            }
    }

    /* initialize */
    tos = bos = NULL;
    forAllNodes(i) {
        i->d = WHITE;
        /*   buckets[i-nodes].firstActive = NULL; */
        buckets[i - nodes].firstActive = sentinelNode;
        i->current = i->first;
    }

    /* eliminate flow cycles, topologicaly order vertices */
    forAllNodes(i)if ((i->d == WHITE) && (i->excess > 0) &&
                      (i != source) && (i != sink)) {
            r = i;
            r->d = GREY;
            do {
                for (; i->current != (i + 1)->first; i->current++) {
                    a = i->current;
                    if ((cap[a - arcs] == 0) && (a->resCap > 0)) {
                        j = a->head;
                        if (j->d == WHITE) {
                            /* start scanning j */
                            j->d = GREY;
                            buckets[j - nodes].firstActive = i;
                            i = j;
                            break;
                        } else if (j->d == GREY) {
                            /* find minimum flow on the cycle */
                            delta = a->resCap;
                            while (1) {
                                delta = min (delta, j->current->resCap);
                                if (j == i)
                                    break;
                                else
                                    j = j->current->head;
                            }

                            /* remove delta flow units */
                            j = i;
                            while (1) {
                                a = j->current;
                                a->resCap -= delta;
                                a->rev->resCap += delta;
                                j = a->head;
                                if (j == i)
                                    break;
                            }

                            /* backup DFS to the first saturated arc */
                            restart = i;
                            for (j = i->current->head; j != i; j = a->head) {
                                a = j->current;
                                if ((j->d == WHITE) || (a->resCap == 0)) {
                                    j->current->head->d = WHITE;
                                    if (j->d != WHITE)
                                        restart = j;
                                }
                            }

                            if (restart != i) {
                                i = restart;
                                i->current++;
                                break;
                            }
                        }
                    }
                }

                if (i->current == (i + 1)->first) {
                    /* scan of i complete */
                    i->d = BLACK;
                    if (i != source) {
                        if (bos == NULL) {
                            bos = i;
                            tos = i;
                        } else {
                            i->bNext = tos;
                            tos = i;
                        }
                    }

                    if (i != r) {
                        i = buckets[i - nodes].firstActive;
                        i->current++;
                    } else
                        break;
                }
            } while (1);
        }


    /* return excesses */
    /* note that sink is not on the stack */
    if (bos != NULL) {
        for (i = tos; i != bos; i = i->bNext) {
            a = i->first;
            while (i->excess > 0) {
                if ((cap[a - arcs] == 0) && (a->resCap > 0)) {
                    if (a->resCap < i->excess)
                        delta = a->resCap;
                    else
                        delta = i->excess;
                    a->resCap -= delta;
                    a->rev->resCap += delta;
                    i->excess -= delta;
                    a->head->excess += delta;
                }
                a++;
            }
        }
        /* now do the bottom */
        i = bos;
        a = i->first;
        while (i->excess > 0) {
            if ((cap[a - arcs] == 0) && (a->resCap > 0)) {
                if (a->resCap < i->excess)
                    delta = a->resCap;
                else
                    delta = i->excess;
                a->resCap -= delta;
                a->rev->resCap += delta;
                i->excess -= delta;
                a->head->excess += delta;
            }
            a++;
        }
    }
}


/* gap relabeling */

int gap(emptyB)
        bucket *emptyB;

{

    bucket *l;
    node *i;
    long r;           /* index of the bucket before l  */
    int cc;          /* cc = 1 if no nodes with positive excess before
		      the gap */

    gapCnt++;
    r = (emptyB - buckets) - 1;

    /* set labels of nodes beyond the gap to "infinity" */
    for (l = emptyB + 1; l <= buckets + dMax; l++) {
        /* this does nothing for high level selection
        for (i = l -> firstActive; i != sentinelNode; i = i -> bNext) {
          i -> d = n;
          gNodeCnt++;
        }
        l -> firstActive = sentinelNode;
        */

        for (i = l->firstInactive; i != sentinelNode; i = i->bNext) {
            i->d = n;
            gNodeCnt++;
        }

        l->firstInactive = sentinelNode;
    }

    cc = (aMin > r) ? 1 : 0;

    dMax = r;
    aMax = r;

    return (cc);

}

/*--- relabelling node i */

long relabel(i)

        node *i;   /* node to relabel */

{

    node *j;
    long minD;     /* minimum d of a node reachable from i */
    arc *minA;    /* an arc which leads to the node with minimal d */
    arc *a;

    assert(i->excess > 0);

    relabelCnt++;
    workSinceUpdate += BETA;

    i->d = minD = n;
    minA = NULL;

    /* find the minimum */
    forAllArcs(i, a) {
        workSinceUpdate++;
        if (a->resCap > 0) {
            j = a->head;
            if (j->d < minD) {
                minD = j->d;
                minA = a;
            }
        }
    }

    minD++;

    if (minD < n) {

        i->d = minD;
        i->current = minA;

        if (dMax < minD) dMax = minD;

    } /* end of minD < n */

    return (minD);

} /* end of relabel */


/* discharge: push flow out of i until i becomes inactive */

void discharge(i)

        node *i;

{

    node *j;                 /* successor of i */
    long jD;                 /* d of the next bucket */
    bucket *lj;               /* j's bucket */
    bucket *l;                /* i's bucket */
    arc *a;                 /* current arc (i,j) */
    cType delta;
    arc *stopA;

    assert(i->excess > 0);
    assert(i != sink);
    do {

        jD = i->d - 1;
        l = buckets + i->d;

        /* scanning arcs outgoing from  i  */
        for (a = i->current, stopA = (i + 1)->first; a != stopA; a++) {
            if (a->resCap > 0) {
                j = a->head;

                if (j->d == jD) {
                    pushCnt++;
                    if (a->resCap < i->excess)
                        delta = a->resCap;
                    else
                        delta = i->excess;
                    a->resCap -= delta;
                    a->rev->resCap += delta;

                    if (j != sink) {

                        lj = buckets + jD;

                        if (j->excess == 0) {
                            /* remove j from inactive list */
                            iDelete(lj, j);
                            /* add j to active list */
                            aAdd(lj, j);
                        }
                    }

                    j->excess += delta;
                    i->excess -= delta;

                    if (i->excess == 0) break;

                } /* j belongs to the next bucket */
            } /* a  is not saturated */
        } /* end of scanning arcs from  i */

        if (a == stopA) {
            /* i must be relabeled */
            relabel(i);

            if (i->d == n) break;
            if ((l->firstActive == sentinelNode) &&
                (l->firstInactive == sentinelNode)
                    )
                gap(l);

            if (i->d == n) break;
        } else {
            /* i no longer active */
            i->current = a;
            /* put i on inactive list */
            iAdd(l, i);
            break;
        }
    } while (1);
}


/* go from higher to lower buckets, push flow */
void wave() {

    node *i;
    bucket *l;

    for (l = buckets + aMax; l > buckets; l--) {
        for (i = l->firstActive; i != sentinelNode; i = l->firstActive) {
            aRemove(l, i);

            assert(i->excess > 0);
            discharge(i);

        }
    }
}


/* first stage  -- maximum preflow*/

void stageOne() {

    node *i;
    bucket *l;             /* current bucket */


#if defined(INIT_UPDATE) || defined(OLD_INIT) || defined(WAVE_INIT)
    globalUpdate ();
#endif

    workSinceUpdate = 0;

#ifdef WAVE_INIT
    wave();
#endif

    /* main loop */
    while (aMax >= aMin) {
        l = buckets + aMax;
        i = l->firstActive;

        if (i == sentinelNode)
            aMax--;
        else {
            aRemove(l, i);

            assert(i->excess > 0);
            discharge(i);

            if (aMax < aMin)
                break;

            /* is it time for global update? */
            if (workSinceUpdate * globUpdtFreq > nm) {
                globalUpdate();
                workSinceUpdate = 0;
            }

        }

    } /* end of the main loop */

    flow = sink->excess;

}


node *decomposePathInternal(node *n, long int *minCap);

void bfs(){
    node *queue[n+3];
    arc *a;

    qInit(n+3);
    source->d=0;
    qEnqueue(source);
    while(!qEmpty){
        node **current;
        qDequeue(current);

        // if(nNode(*current) == nNode(sink)) break;

        forAllArcs(*current, a){
            if ((a->cap == 0) || (a->cap == a->resCap)) continue;
            if (a->head->d == -1){
                if (a->head == sink) {
                    fprintf(stderr, "Linked from %ld with cap=%ld and resCap=%ld\n", nNode(*current), a->cap, a->resCap);
                    fflush(stderr);
                }
                a->head->d = (*current)->d + 1;
                fprintf(stderr, "Discovered node %ld at distance %ld cap=%ld resCap=%ld\n", nNode(a->head), (*current)->d + 1, a->cap, a->resCap);
                qEnqueue(a->head);
            }
        }
    }
    node *u;
    forAllNodes(u) {
        fprintf(stderr, "Node %ld at distance %ld\n", nNode(u), u->d);
    }
}



void hipr(
        long ninput,
        long minput,
        long *tails,
        long *heads,
        long *weights,
        long s,
        long t,
        int **output_set,
        long **mheads,
        long **mtails,
        long **mweights,
        long *nedges,
        long *fflow,
        int route_flag)
{
#if (defined(PRINT_FLOW) || defined(CHECK_SOLUTION))
    node *i;
    arc *a;
#endif

#ifdef PRINT_FLOW
    long ni, na;
#endif
    node *j;
    int cc;
#ifdef CHECK_SOLUTION
    excessType sum;
    bucket *l;
#endif
    /* fprintf(stderr,"calling hi_pr\n"); */


    /* t1 = timer();

   t2 = t1;
    */


    loadflowproblem(ninput, minput, tails, heads, weights, s, t, &n, &m, &nodes, &arcs, &cap, &source, &sink, &nMin);


    /*  fprintf(stderr,"nodes:       %10ld\nc arcs:        %10ld\nc\n", n, m); */

    cc = allocDS();
    if (cc) {
        fprintf(stderr, "Allocation error\n");
        exit(1);
    }

    init();
    stageOne();

    /*
    t2 = timer() - t2;
    */

    /*  fprintf (stderr,"c flow:       %12.01f\n", flow); */

#ifndef CUT_ONLY
    stageTwo();

    /*
    t1 = timer() - t1;

    fprintf (stderr,"\nc time:        %10.2f\n", t1);
    */
#endif
    /*
    fprintf (stderr,"c cut tm:      %10.2f\n", t2);
    */


#ifdef CHECK_SOLUTION

    /* check if you have a flow (pseudoflow) */
    /* check arc flows */
    forAllNodes(i) {
        forAllArcs(i, a) {
            if (cap[a - arcs] > 0) /* original arc */
                if ((a->resCap + a->rev->resCap != cap[a - arcs])
                    || (a->resCap < 0)
                    || (a->rev->resCap < 0)) {
                    printf("ERROR: bad arc flow\n");
                    exit(2);
                }
        }
    }

    /* check conservation */
    forAllNodes(i)if ((i != source) && (i != sink)) {
#ifdef CUT_ONLY
            if (i->excess < 0) {
          printf("ERROR: nonzero node excess\n");
          exit(2);
            }
#else
            if (i->excess != 0) {
                printf("ERROR: nonzero node excess\n");
                exit(2);
            }
#endif

            sum = 0;
            forAllArcs(i, a) {
                if (cap[a - arcs] > 0) /* original arc */
                    sum -= cap[a - arcs] - a->resCap;
                else
                    sum += a->resCap;
            }

            if (i->excess != sum) {
                printf("ERROR: conservation constraint violated\n");
                exit(2);
            }
            }

    /* check if mincut is saturated */
    aMax = dMax = 0;
    for (l = buckets; l < buckets + n; l++) {
        l->firstActive = sentinelNode;
        l->firstInactive = sentinelNode;
    }
    globalUpdate();
    if (source->d < n) {
        printf("ERROR: the solution is not optimal\n");
        exit(2);
    }

    /*  fprintf(stderr,"c\nc Solution checks (feasible and optimal)\nc\n");*/
#endif

#ifdef PRINT_STAT
    printf ("c pushes:      %10ld\n", pushCnt);
    printf ("c relabels:    %10ld\n", relabelCnt);
    printf ("c updates:     %10ld\n", updateCnt);
    printf ("c gaps:        %10ld\n", gapCnt);
    printf ("c gap nodes:   %10ld\n", gNodeCnt);
    printf ("c\n");
#endif

#ifdef PRINT_FLOW
    printf ("c flow values\n");
    forAllNodes(i) {
      ni = nNode(i);
      forAllArcs(i,a) {
    na = nArc(a);
    if ( cap[na] > 0 )
      printf ( "f %7ld %7ld %12ld\n",
          ni, nNode( a -> head ), cap[na] - ( a -> resCap )
          );
      }
    }
    printf("c\n");
#endif

#ifdef PRINT_CUT
    globalUpdate();
    printf ("c nodes on the sink side\n");
    forAllNodes(j)
      if (j->d < n)
        printf("c %ld\n", nNode(j));

#endif
    /***************************************************************
     BEGINNING OF MODIFICATIONS
     ***************************************************************
     ***************************************************************/

    /*  t1 = timer();*/
    /* RETRIEVE CUT FOUND - CODE BY REID */

    *output_set = calloc(sizeof(**output_set), n - 2);

    forAllNodes(j)if (j->d < n) {
            if (nNode(j) - nMin < n - 2) {
                (*output_set)[nNode(j) - nMin] = 1;
            }
        }

    /* RETRIEVE FLOW */
    *fflow = (long) flow;

    /* RETRIEVE ROUTED GRAPH - CODE BY SATISH */

    if (route_flag == 1) {
        fprintf(stderr, "Ready for matching\n");
        fflush(stderr);
        long int matchingCapacity;
        long *reallocPtr;

        int k;

        source->excess = 0;
        forAllArcs(source, a) {
            long na = nArc(a);
            if (cap[na] == 0) continue;
            source->excess += cap[na] - a->resCap;
        }
        /* INITIALIZE MATCHING ARRAYS */

        matchingCapacity = 0;
        k = 0;

        node* mhead;
        node* mtail;
        long mweight;
        dynamic_tree_t *p;
        while(source->excess != 0) {
            fprintf(stderr, "source->excess=%llu\n", source->excess);
            fflush(stderr);
            long sink_excess = 0;
            forAllNodes(i) {
                forAllArcs(i, a) {
                    long na = nArc(a);
                    if ((a->cap > 0) && (a->head == sink)) {
                        sink_excess += cap[na] - a->resCap;
                    }
                }
            }
            fprintf(stderr, "sink->excess=%ld\n", sink_excess);
            fflush(stderr);

            sleep(2);

            forAllNodes(i) {
                i->d = -1;
                i->current = i->first;
            }
            bfs();

            p = dt_init(n, nodes, source);
            while (source->current < (source + 1)->first) {
                while (p->cur_node != sink) {
                    int link_flag = 0;
                    for (; p->cur_node->current < (p->cur_node + 1)->first; p->cur_node->current++) {   // Find suitable edge or exhaust edges
                        arc* cur_arc = p->cur_node->current;
                        if (cap[nArc(cur_arc)] == 0) continue;              // Reverse arc, not important.

                        if (cur_arc->head == sink) {
                            fprintf(stderr, "Found edge %ld -> %ld. cap=%ld. resCap=%ld. d=%ld. next d=%ld\n",
                                    nNode(p->cur_node), nNode(cur_arc->head), cur_arc->cap, cur_arc->resCap, p->cur_node->d, cur_arc->head->d);
                            fflush(stderr);
                            getchar();
                        }
                        if ((cur_arc->cap == cur_arc->resCap) ||
                            (p->cur_node->d + 1 != cur_arc->head->d)) {
                            continue;
                        }
                        p->cur_node->current++; // added by Xifan
                        fprintf(stderr, "Linking %ld -> %ld\n", nNode(p->cur_node), nNode(cur_arc->head));
                        fflush(stderr);
                        dt_link(p, p->cur_node, cur_arc->head, cur_arc);
                        link_flag = 1;
                        break;
                    }
                    if (link_flag == 1) {
                        continue;
                    }
                    if (p->cur_node->current == (p->cur_node + 1)->first) {             // if no suitable edges cut tail
                        node* previous;
                        if ((previous = dt_before(p, p->cur_node)) != NULL) {              // Checks that a previous exists.
                                                                                        // The alternative is that p is the source
                            fprintf(stderr, "About to cut to %ld\n", nNode(previous));
                            fflush(stderr);
                            dt_cut(p, previous);
                            // following line commented out by Xifan
                            // p->cur_node->current++; // TODO: this update is problematic. It should be previous->current++
                        } else {
                            break;
                        }
                    }
                }

                if (p->cur_node != sink) {
                    break;
                }

                dt_findPath(p, &mhead, &mtail, &mweight);

                fprintf(stderr, "Current mweight = %ld\n", mweight);
                fflush(stderr);

                if(!mweight) continue;

                source->excess -= mweight;

                if (k >= matchingCapacity) {
                    if (!matchingCapacity) matchingCapacity = 2 * n;
                    else matchingCapacity = 2 * matchingCapacity;
                    reallocPtr = *mheads;
                    *mheads = realloc(*mheads, sizeof(**mheads) * matchingCapacity);
                    if (NULL == *mheads) {
                        free(reallocPtr);
                        fprintf(stderr, "Failed to allocate mheads for %ld places\n", matchingCapacity);
                        exit(1);
                    }

                    reallocPtr = *mtails;
                    *mtails = realloc(*mtails, sizeof(**mtails) * matchingCapacity);
                    if (NULL == *mtails) {
                        free(reallocPtr);
                        fprintf(stderr, "Failed to allocate mheads for %ld places\n", matchingCapacity);
                        exit(1);
                    }

                    reallocPtr = *mweights;
                    *mweights = realloc(*mweights, sizeof(**mweights) * matchingCapacity);
                    if (NULL == *mweights) {
                        free(reallocPtr);
                        fprintf(stderr, "Failed to allocate mheads for %ld places\n", matchingCapacity);
                        exit(1);
                    }
                }

                (*mheads)[k] = nNode(mhead);
                (*mtails)[k] = nNode(mtail);
                (*mweights)[k] = mweight;

                (*mtails)[k + 1] = nNode(mhead);
                (*mheads)[k + 1] = nNode(mtail);
                (*mweights)[k + 1] = mweight;

                k = k + 2;

            }
            dt_cleanUp(p);
        }
        *nedges = k;
    }


    /*  fprintf(stderr, "rem tm: %f//\n",  timer() - t1);		*/


    /* Free data structures */
    free(nodes - nMin);               /* address of the array of nodes */ /*MEMORY LEAK*/
    free(arcs);             /* address of the array of arcs */
    free(cap);              /* address of the array of capasities */
    free(buckets);              /* address of the array of capasities */
}

int loadflowproblem(n, m, tails, heads, weights, s, t,
                    n_ad, m_ad, nodes_ad, arcs_ad, cap_ad,
                    source_ad, sink_ad, node_min_ad)
/* input */
        long n;
        long m;
        long *tails;
        long *heads;
        long *weights;
        long s;
        long t;

/* output */
        long *n_ad;                 /* address of the number of nodes */
        long *m_ad;                 /* address of the number of arcs */
        node **nodes_ad;            /* address of the array of nodes */
        arc **arcs_ad;             /* address of the array of arcs */
        long **cap_ad;              /* address of the array of capasities */
        node **source_ad;           /* address of the pointer to the source */
        node **sink_ad;             /* address of the pointer to the source */
        long *node_min_ad;          /* address of the minimal node */
{

    long
            node_min = 0,             /* minimal no of node  */
    node_max = 0,             /* maximal no of node */
    *arc_first = NULL,         /* internal array for holding
                                     - node degree
                                     - position of the first outgoing arc */
    *arc_tail = NULL,          /* internal array: tails of the arcs */
    source = 0,               /* no of the source */
    sink = 0,                 /* no of the sink */
    /* temporary variables carrying no of nodes */
    head, tail, i;

    long
    /* temporary variables carrying no of arcs */
    last, arc_num, arc_new_num;

    node *nodes = NULL,            /* pointer to the node structure */
    *head_p,
            *ndp;

    arc *arcs = NULL,             /* pointer to the arc structure */
    *arc_current = NULL,
            *arc_new,
            *arc_tmp;

    long *acap = NULL,             /* array of capasities */
    cap;                    /* capasity of the current arc */

    long
            pos_current = 0;          /* 2*no_alines */

    long k;                      /* temporary */

/* allocating memory for  'nodes', 'arcs'  and internal arrays */
    nodes = (node *) calloc(n + 2, sizeof(node));
    arcs = (arc *) calloc(2 * m + 1, sizeof(arc));
    arc_tail = (long *) calloc(2 * m, sizeof(long));
    arc_first = (long *) calloc(n + 2, sizeof(long));
    acap = (long *) calloc(2 * m, sizeof(long));
    /* arc_first [ 0 .. n+1 ] = 0 - initialized by calloc */

    /* setting pointer to the first arc */
    arc_current = arcs;

    source = s;
    sink = t;

    node_max = 0;
    node_min = n;

    for (k = 0; k < m; k++) {
        tail = tails[k];
        head = heads[k];
        cap = weights[k];

        /* fprintf(stderr,"%d of %d: %d %d %d\n",k,m,tail,head,cap); */

        /* no of arcs incident to node i is stored in arc_first[i+1] */
        arc_first[tail + 1]++;
        arc_first[head + 1]++;

        /* storing information about the arc */
        arc_tail[pos_current] = tail;
        arc_tail[pos_current + 1] = head;
        arc_current->head = nodes + head;
        arc_current->resCap = cap;
        arc_current->cap = cap;
        arc_current->rev = arc_current + 1;
        (arc_current + 1)->head = nodes + tail;
        (arc_current + 1)->resCap = 0;
        (arc_current + 1)->cap = 0;
        (arc_current + 1)->rev = arc_current;


        /* searching minimum and maximum node */
        if (head < node_min) node_min = head;
        if (tail < node_min) node_min = tail;
        if (head > node_max) node_max = head;
        if (tail > node_max) node_max = tail;

        arc_current += 2;
        pos_current += 2;
    }

/********** ordering arcs - linear time algorithm ***********/

/* first arc from the first node */
    (nodes + node_min)->first = arcs;

/* before below loop arc_first[i+1] is the number of arcs outgoing from i;
   after this loop arc_first[i] is the position of the first
   outgoing from node i arcs after they would be ordered;
   this value is transformed to pointer and written to node.first[i]
   */

    for (i = node_min + 1; i <= node_max + 1; i++) {
        arc_first[i] += arc_first[i - 1];
        (nodes + i)->first = arcs + arc_first[i];
    }


    for (i = node_min; i < node_max; i++) /* scanning all the nodes
                                            exept the last*/
    {

        last = ((nodes + i + 1)->first) - arcs;
        /* arcs outgoing from i must be cited
         from position arc_first[i] to the position
         equal to initial value of arc_first[i+1]-1  */

        for (arc_num = arc_first[i]; arc_num < last; arc_num++) {
            tail = arc_tail[arc_num];

            while (tail != i)
                /* the arc no  arc_num  is not in place because arc cited here
                   must go out from i;
                   we'll put it to its place and continue this process
                   until an arc in this position would go out from i */

            {
                arc_new_num = arc_first[tail];
                arc_current = arcs + arc_num;
                arc_new = arcs + arc_new_num;

                /* arc_current must be cited in the position arc_new
                   swapping these arcs:                                 */

                head_p = arc_new->head;
                arc_new->head = arc_current->head;
                arc_current->head = head_p;

                cap = arc_new->resCap;
                arc_new->resCap = arc_current->resCap;
                arc_current->resCap = cap;

                cap = arc_new->cap;
                arc_new->cap = arc_current->cap;
                arc_current->cap = cap;

                if (arc_new != arc_current->rev) {
                    arc_tmp = arc_new->rev;
                    arc_new->rev = arc_current->rev;
                    arc_current->rev = arc_tmp;

                    (arc_current->rev)->rev = arc_current;
                    (arc_new->rev)->rev = arc_new;
                }

                arc_tail[arc_num] = arc_tail[arc_new_num];
                arc_tail[arc_new_num] = tail;

                /* we increase arc_first[tail]  */
                arc_first[tail]++;

                tail = arc_tail[arc_num];
            }
        }
        /* all arcs outgoing from  i  are in place */
    }

/* -----------------------  arcs are ordered  ------------------------- */

/*----------- constructing lists ---------------*/


    for (ndp = nodes + node_min; ndp <= nodes + node_max; ndp++)
        ndp->first = (arc *) NULL;

    for (arc_current = arcs + (2 * m - 1); arc_current >= arcs; arc_current--) {
        arc_num = arc_current - arcs;
        tail = arc_tail[arc_num];
        ndp = nodes + tail;
        /* avg
        arc_current -> next = ndp -> first;
        */
        ndp->first = arc_current;
    }


/* ----------- assigning output values ------------*/
    *m_ad = m;
/*n_ad = node_max - node_min + 1;*/
    *n_ad = n;
    *source_ad = nodes + source;
    *sink_ad = nodes + sink;
    *node_min_ad = node_min;
    *nodes_ad = nodes + node_min;
    *arcs_ad = arcs;
    *cap_ad = acap;

    for (arc_current = arcs, arc_num = 0;
         arc_num < 2 * m;
         arc_current++, arc_num++
            )
        acap[arc_num] = arc_current->resCap;

    if (source < node_min || source > node_max)
        fprintf(stderr, "bad value for source");

    if ((*source_ad)->first == (arc *) NULL ||
        (*sink_ad)->first == (arc *) NULL)
        fprintf(stderr, "no arc goes out of the source");

/* free internal memory */
    free(arc_first);
    free(arc_tail);

/* Thanks God! all is done */
    return (0);

}


node *decomposePathInternal(node *n, long int *minCap) {
    node *i;
    arc *a;
    arc *stopA;

    if (n->d < 0) {
        if (n->d == -2) {
            printf("Cycle detection software error.\n");
        }
        n->d = -2;
        return NULL;
    }

    n->d = -1; /* node on stack. node. */

    if (n == sink) {
        printf("Returning sink, means direct source-sink path.\n");
        n->d = 0;
        return n;
    }

    forAllArcs(n, a) {
        int na = nArc(a);
        if (cap[na] > 0) {
            while (a->resCap < cap[na]) {
                int thisCap = cap[na] - a->resCap;
                int tempRes = a->resCap;

                if (*minCap > thisCap)
                    *minCap = thisCap;
                if (a->head == sink) {
                    a->resCap += *minCap;
                    n->d = 0;
                    return n;
                }

                i = decomposePathInternal(a->head, minCap);
                a->resCap = tempRes + *minCap;
                if (((i == NULL) && (n->d == -2))) {
                    if (i == NULL)
                        n->d = -1;
                } else {
                    /*
                    if (i==NULL)
                      {
                        printf ("Cancelling cycle\n");
                      }
                    */
                    n->d = 0;
                    return i;
                }
            }
        }
    }
    printf("ERROR: should never call decomposePathInternal on node with no outgoing flow!\n");
    n->d = 0;
    return n;
}
