/* 
 * GraphAlpha.h - GraphAlpha v1.3
 * =============================================================
 * Header file of GraphAlpha (library of graph algorithms)
 */

#ifndef GraphAlpha_H
#define GraphAlpha_H

/************************** DATA STRUCTURES *************************/
typedef struct c_arc {                  /* data type of arc */
	unsigned int index;									  /* index of neighbor */
	double weight;											  /* weight of respective arc */
  struct c_arc* next;									  /* next neighbor of the source node */
} arc;

typedef struct c_node {                 /* data type of node */
	arc* successors;										  /* linked list of successors */
} node;

typedef struct c_graph {                /* data type of graph */
  unsigned int size;                    /* # of nodes */
  unsigned int order;                   /* # of arcs */

	node* nodes;												  /* array of nodes of graph */
} graph;

typedef struct c_heap {                 /* data type of k-regular heap */
  unsigned int regularity;              /* regularity of the heap */
  int* indices;                         /* array of index-mapping from vertices
                                           to heap elements */
  int* heap_arr;                        /* array of the heap */
  int bottom;                           /* end of heap in heap_arr */
} heap;

typedef struct {                        /* data type of a single bit */
  unsigned int bit : 1;
} c_bit;

typedef struct c_llist_node {           /* data type of a node of queue */
  unsigned int val;
  struct c_llist_node* next;
} llist_node;

typedef struct c_queue {                /* data type of the whole queue */
  llist_node* head;
  llist_node* tail;
} queue;
/************************** DATA STRUCTURES *************************/

/************************** FUNCTION HEADERS ************************/
double** matrixConvert(graph* p_gr);    /* receives pointer to list of
                                           successor, gives out adjacancy
                                           matrix
                                           */
double** symMtrx(graph* p_gr);          /* symetrization of arcs, returns a
                                           matrix */
double* Dijkstra (int v, graph* p_gr);  /* Dijkstra's algorithm for shortest
                                           paths from initial point */
/* receives adjacency matrix (with distances), gives out list of successors */
graph* listConvert(double** matrix, int dim);
graph* listConvertBit(c_bit** matrix, int dim);
double** FloydWarshall(graph* p_gr);    /* Floyd-Warshall's algorithm for
                                           finding shortest paths between every
                                           2 nodes */
void Tarjan(graph* p_gr);               /* Tarjan's algorithm for strongly
                                           connected components */
int TSort(graph* p_gr, int* seq);       /* topological ordering */
c_bit** transitiveClosure(graph* p_gr); /* transitive closure */
c_bit** transitiveReduction(graph* p_gr);/* transitive reduction */
/************************** FUNCTION HEADERS ************************/

#endif
