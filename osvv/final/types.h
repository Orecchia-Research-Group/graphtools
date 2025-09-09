/* defs.h */

#ifndef FILE_TYPES_SEEN
#define FILE_TYPES_SEEN

#include <stdint.h>

typedef uint64_t excessType;

/*typedef unsigned long cType;*/
typedef int64_t cType;

typedef  /* arc */
   struct arcSt
{
   cType           cap;             /* maximum capacity */
   cType           resCap;          /* residual capacity */
   struct nodeSt   *head;           /* arc head */
   struct arcSt    *rev;            /* reverse arc */
} arc;

typedef  /* node */
   struct nodeSt
{
   arc             *first;           /* first outgoing arc */
   arc             *current;         /* current outgoing arc */
   excessType      excess;           /* excess at the node 
				        change to double if needed */
   int64_t         d;                /* distance label */
   struct nodeSt   *bNext;           /* next node in bucket */
   struct nodeSt   *bPrev;           /* previous node in bucket */
} node;


typedef /* bucket */
   struct bucketSt
{
  node             *firstActive;      /* first node with positive excess */
  node             *firstInactive;    /* first node with zero excess */
} bucket;

#endif
