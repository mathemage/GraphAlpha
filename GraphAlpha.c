/* 
 * GraphAlpha.h - GraphAlpha v1.3
 * =============================================================
 * Functions of GraphAlpha (library of graph algorithms)
 */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <malloc.h>
#include "GraphAlpha.h"

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  maxAbs
 *  Description:  returns parameter with greater absolute value
 * =====================================================================================
 */
double maxAbs(double a, double b)
{
  return (a*a > b*b) ? a : b;
}		/* -----  end of function maxAbs  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  symMtrx
 *  Description:  returns matrix of graph with symmetrized arcs; when 2
 *  non-null weights are present, the one with greater absolute value is chosen
 * =====================================================================================
 */
double ** symMtrx (graph* p_gr)
{
  double** matrix = matrixConvert(p_gr);
  int i, j;
  int res_value;
  for ( i = 0; i < p_gr->size; i += 1 ) {
    for ( j = i; j < p_gr->size; j += 1 ) {
      if (DBL_MAX == matrix[i][j])
        matrix[i][j] = matrix[j][i];
      else if (DBL_MAX == matrix[j][i])
        matrix[j][i] = matrix[i][j];
      else
        matrix[i][j] = matrix[j][i] = maxAbs(matrix[i][j], matrix[j][i]);
    }
  }
  return matrix;
}		/* -----  end of function symMtrx  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  swap
 *  Description:  swap two int values
 * =====================================================================================
 */
void swap(int* arr, int idx1, int idx2)
{
  int tmp = arr[idx1];
  arr[idx1] = arr[idx2];
  arr[idx2] = tmp;
}		/* -----  end of function swap  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  heapDecrease
 *  Description:  
 * =====================================================================================
 */
void heapDecrease(heap* p_h, double* distances, int v, double value)
{
  /* graph-index of v's father */
  int idx_dad = v;                                               /* root has no daddy */
  if ( p_h->indices[v] > 0 ) {                                   /* otherwise */
    idx_dad = p_h->heap_arr[ (p_h->indices[v]-1)/p_h->regularity ];  
  }

  distances[v] = value;

  /* rotate upwards */
  while (p_h->indices[idx_dad] < p_h->indices[v]
      && distances[idx_dad] > distances[v]) {
    swap(p_h->heap_arr, p_h->indices[idx_dad], p_h->indices[v]); /* swap in heap */
    swap(p_h->indices, idx_dad, v);                              /* swap of indices */

    v = idx_dad;

    /* parent of our parent */
    idx_dad = p_h->heap_arr[v];                                  /* root has no daddy */
    if ( p_h->indices[v] > 0 ) {                                 /* otherwise */
      idx_dad = p_h->heap_arr[ (p_h->indices[v]-1)/p_h->regularity ];  
    }
  }
}		/* -----  end of function heapDecrease  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  heapInsert
 *  Description:  insert new element into heap
 * =====================================================================================
 */
void heapInsert(heap* p_h, double* distances, int v, double value)
{
  p_h->heap_arr[++p_h->bottom] = v;
  p_h->indices[v] = p_h->bottom;
  heapDecrease(p_h, distances, v, value);
}		/* -----  end of function heapInsert  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  heapExtractMin
 *  Description:  remove the root (i.e., smallest element) from the heap
 * =====================================================================================
 */
int heapExtractMin(heap* p_h, double* distances)
{
  int min = p_h->heap_arr[0];

  /* REMOVE THE TOP OF HEAP BY REPLACING WITH THE BOTTOM */
  swap(p_h->indices, p_h->heap_arr[0], p_h->heap_arr[p_h->bottom]);
  swap(p_h->heap_arr, 0, p_h->bottom);
  p_h->bottom--;

  /* ROTATE DOWNWARDS THE PSEUDO-ROOT */
  int i;

  /* heap-indices of father & son */
  int dad = 0;
  /* finding the smallest son */
  int son = 1;
  for ( i = 2; (i <= p_h->regularity)
      && (dad * p_h->regularity + i < p_h->bottom); i += 1 ) {
    if (distances[p_h->heap_arr[son]] > distances[p_h->heap_arr[dad * p_h->regularity + i]]) {
      son = dad * p_h->regularity + i;
    }
  }

  /* rotate downwards */
  while (son <= p_h->bottom
      && distances[p_h->heap_arr[dad]] > distances[p_h->heap_arr[son]]) {
    swap(p_h->indices, p_h->heap_arr[dad], p_h->heap_arr[son]);  /* swap of indices */
    swap(p_h->heap_arr, dad, son);                               /* swap in heap */

    /* new heap-indices of father & son */
    int dad = son;
    /* finding the smallest son */
    int son = dad * p_h->regularity + 1;
    for ( i = 2; (i <= p_h->regularity)
        && (dad * p_h->regularity + i < p_h->bottom); i += 1 ) {
      if (distances[p_h->heap_arr[son]] > distances[p_h->heap_arr[dad * p_h->regularity + i]]) {
        son = dad * p_h->regularity + i;
      }
    }
  }

  return min;
}		/* -----  end of function heapExtractMin  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Dijkstra
 *  Description:  Dijkstra's algorithm for shortest paths from initial point
 * =====================================================================================
 */
double* Dijkstra (int v, graph* p_gr)
{
  /* INIT */
  int i;

  /* the distance-array */
  double* distances = (double*) calloc(p_gr->size, sizeof(double));
  for ( i = 0; i < p_gr->size; i += 1 ) {
    distances[i] = DBL_MAX;
  }
  distances[v] = 0;

  /* the heap */
  heap h;
  h.bottom = -1;
  h.heap_arr = (int*) calloc(p_gr->size, sizeof(int));
  h.indices = (int*) calloc(p_gr->size, sizeof(int));
  for ( i = 0; i < p_gr->size; i += 1 ) {
    h.indices[i] = -1;
  }
  h.regularity = p_gr->order / p_gr->size;
  if (h.regularity < 2) {
    h.regularity = 2;
  }

  heapInsert(&h, distances, v, 0);        /* insert initial graph-node */

  /* DIJKSTRA'S ALGORITHM */
  while (h.bottom >= 0) {                 /* while heap is non-empty */
    v = heapExtractMin(&h, distances);
    
    /* relaxing the arcs leading to neighbors */
    arc* p_cur_arc = p_gr->nodes[v].successors;
    while (NULL != p_cur_arc) {
      /* checking negative weights of edges */
      if (p_cur_arc->weight < 0 ) {
        printf ("Negative edges! Results ungranted...\n");
      }

      /* relax! */
      if (distances[p_cur_arc->index] > distances[v] + p_cur_arc->weight) {
        if (DBL_MAX == distances[p_cur_arc->index]) {    /* if node not visited yet */
          /* insert it */
          heapInsert(&h, distances, p_cur_arc->index, distances[v] + p_cur_arc->weight);
        }
        else {
          /* otherwise directly just decrease it */
          heapDecrease(&h, distances, p_cur_arc->index, distances[v] + p_cur_arc->weight);
        }
      }

      p_cur_arc = p_cur_arc->next;
    }
  }

  return distances;
}		/* -----  end of function Dijkstra  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  double** FloydWarshall(graph* p_gr)
 *  Description:  Floyd-Warshall's algorithm for finding shortest paths between
 *  every 2 nodes
 * =====================================================================================
 */
double** FloydWarshall(graph* p_gr) {
  /* INIT OF MATRIX WITH 0'S ON THE DIAGONAL & INFINITIES ANYWHERE ELSE */
  int i, j, k;
  double** D = (double **) calloc(p_gr->size, sizeof(double *));
  for (i = 0; i < p_gr->size; i++) {
    D[i] = (double *) calloc(p_gr->size, sizeof(double));
    for (j = 0; j < p_gr->size; j++)
      D[i][j] = (i == j) ? 0 : DBL_MAX;
  }

  /* FILL D WITH ACTUAL WEIGHTS */
  arc* p_cur_arc;
  for (i = 0; i < p_gr->size; i++) {
    p_cur_arc = p_gr->nodes[i].successors;              /* neighbors of i-th node */
    
    while (NULL != p_cur_arc) {
      D[i][p_cur_arc->index] = p_cur_arc->weight;
      p_cur_arc = p_cur_arc->next;
    }
  }

  /* THE FLOYD-WARSHALL */
  for (k = 0; k < p_gr->size; k++)
    for (i = 0; i < p_gr->size; i++)
      for (j = 0; j < p_gr->size; j++)
        if (D[i][j] > D[i][k] + D[k][j])
          D[i][j] = D[i][k] + D[k][j];

  return D;
}		/* -----  end of function double** FloydWarshall(graph* p_gr)  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  push
 *  Description:  enqueue into an int-queue
 * =====================================================================================
 */
void push(int value, queue* q)
{
  if ( NULL == q->tail ) {
    q->head = q->tail = (llist_node *) malloc(sizeof(llist_node));
  }
  else {
    q->tail->next = (llist_node *) malloc(sizeof(llist_node));
    q->tail = q->tail->next;
  }
  q->tail->val = value;
}		/* -----  end of function push  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  pop
 *  Description:  dequeue from an int-queue
 * =====================================================================================
 */
int pop(queue* q)
{
  if ( NULL == q->head ) {
    return -1;
  }
  else {
    int res = q->head->val;

    if (NULL == (q->head = q->head->next))
      q->tail = NULL;
    return res;
  }
}		/* -----  end of function pop  ----- */

llist_node* stack = NULL;
int timer = 0;
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  DFS
 *  Description:  DFS with modification for Tarjan's algorithm for find
 *  strongly connected components
 * =====================================================================================
 */
void DFS(graph* p_gr, int v, int* in, int* low, c_bit* visited, c_bit* stacked)
{
  low[v] = in[v];

  /* calculating low values of successors */
  arc* w = p_gr->nodes[v].successors;
  while (NULL != w) {
    if (0 == visited[w->index].bit) {
      visited[w->index].bit = 1;

      /* pushing into stack */
      llist_node* new_node = (llist_node*) malloc(sizeof(llist_node));
      new_node->val = w->index;
      new_node->next = stack;
      stack = new_node;
      stacked[w->index].bit = 1;

      in[w->index] = ++timer;
      DFS(p_gr, w->index, in, low, visited, stacked);
      if (low[v] > low[w->index]) {                /* improvement via tree arc */
        low[v] = low[w->index];
      }
    }
    if (1 == stacked[w->index].bit && in[w->index] < in[v]
        && low[v] > in[w->index]) {                /* improvement via reverse arc */
      low[v] = in[w->index];
    }
    w = w->next;
  }

  if (low[v] == in[v]) {
    /* popping nodes with in[w] >= in[v] from stack */
    llist_node* del;
    printf("( ");
    while (NULL != stack && in[stack->val] >= in[v]) {
      printf("%d ", stack->val + 1);
      stacked[stack->val].bit = 0;
      del = stack;
      stack = stack->next;
      free(del);
    }
    printf(")\n");
  }
}		/* -----  end of function DFS  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Tarjan
 *  Description:  Tarjan's algorithm for strongly connected components 
 * =====================================================================================
 */
void Tarjan(graph* p_gr)
{
  /* INIT */
  int in[p_gr->size]; int low[p_gr->size]; c_bit stacked[p_gr->size];
  stack = NULL;
  timer = 0;

  int v;
  c_bit visited[p_gr->size];
  for ( v = 0; v < p_gr->size; v += 1 ) {
    visited[v].bit = 0;
    stacked[v].bit = 0;
  }

  /* SOMEWHERE TO BEGIN */
  llist_node* new_node;
  for ( v = 0; v < p_gr->size; v += 1 ) {
    if (0 == visited[v].bit) {
      visited[v].bit = 1;

      new_node = (llist_node *) malloc (sizeof(llist_node));
      new_node->val = v;
      new_node->next = stack;
      stack = new_node;
      stacked[v].bit = 1;

      in[v] = ++timer;
      DFS(p_gr, v, in, low, visited, stacked);
    }
  }
}		/* -----  end of function Tarjan  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  TSort
 *  Description:  topological ordering -- receives graph & empty array, if
 *  possible returns topological ordering of graph in the array otherwise
 *  returns code -1 in case of non-acyclic graphs
 * =====================================================================================
 */
int TSort(graph* p_gr, int* seq)
{
  /* CREATING ARRAY WITH COUNTS OF NODES' PREDECESSORS */
  unsigned int cnt_predec[p_gr->size];

  /* init */
  int i;
  for ( i = 0; i < p_gr->size; i += 1 ) {
    cnt_predec[i] = 0;
  }

  /* fill with the counts */
  arc* p_cur_arc;
  for ( i = 0; i < p_gr->size; i += 1 ) {
    p_cur_arc = p_gr->nodes[i].successors;

    while ( NULL != p_cur_arc ) {
      cnt_predec[p_cur_arc->index]++;
      p_cur_arc = p_cur_arc->next;
    }
  }

  /* TEARING OFF SOURCE-NODES */
  queue* q = (queue *) malloc(sizeof(queue));
  q->head = q->tail = NULL;
  int cnt_torn = 0;             /* count of torn-off nodes */

  /* searching for original source-nodes*/
  for ( i = 0; i < p_gr->size; i += 1 ) {
    if ( 0 == cnt_predec[i] ) {
      push(i, q);
    }
  }

  int idx;
  int j = 0;                                    /* index in seq[] */
  /* seeking another source-nodes */
  while ( -1 != (idx = pop(q)) ) {              /* until empty queue */
    seq[j++] = idx;                             /* insert into top. ordering */

    p_cur_arc = p_gr->nodes[idx].successors;    /* decrease neigbors' count of
                                                   predecessors */
    while ( NULL != p_cur_arc ) {                  
      if (0 == --cnt_predec[p_cur_arc->index])  /* new source node */
        push(p_cur_arc->index, q);
      p_cur_arc = p_cur_arc->next;
    }
  }
  
  if ( j == p_gr->size ) {
    return 0;                                   /* success! */
  }
  else {
    return -1;                                  /* failure */
  }
}		/* -----  end of function TSort  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  transitiveClosure
 *  Description:  transitive closure -- modification of Floyd-Warshall's
 *  algorithm, instead of evaluating distance, only accessibility is taken into
 *  account 
 * =====================================================================================
 */
c_bit** transitiveClosure(graph* p_gr)
{
  /* INIT OF MATRIX WITH 1'S ON THE DIAGONAL & 0'S ANYWHERE ELSE */
  /* nodes are always accessible to themselves */
  int i, j, k;
  c_bit** D = (c_bit **) calloc(p_gr->size, sizeof(c_bit *));
  for (i = 0; i < p_gr->size; i++) {
    D[i] = (c_bit *) calloc(p_gr->size, sizeof(c_bit));
    for (j = 0; j < p_gr->size; j++)
      D[i][j].bit = (i == j) ? 1 : 0;
  }

  /* FILL D WITH ACCESSIBILITIES THROUGH ARCS */
  arc* p_cur_arc;
  for (i = 0; i < p_gr->size; i++) {
    p_cur_arc = p_gr->nodes[i].successors;       /* neighbors of i-th node */

    while (NULL != p_cur_arc) {
      D[i][p_cur_arc->index].bit = 1;
      p_cur_arc = p_cur_arc->next;
    }
  }

  /* THE MODIFIED FLOYD-WARSHALL */
  for (k = 0; k < p_gr->size; k++)
    for (i = 0; i < p_gr->size; i++)
      for (j = 0; j < p_gr->size; j++)
        D[i][j].bit |= D[i][k].bit & D[k][j].bit; /* true if there's already
                                                     i->j path or there are
                                                     both i->k and k->j paths
                                                     */

  return D;
}		/* -----  end of function transitiveClosure  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  transitiveReduction
 *  Description:  transitive reduction -- if possible: R-=R \ (R ° R+) - viz
 *  http://en.wikipedia.org/wiki/Transitive_reduction#Graph_algorithms_for_transitive_reduction
 *  if not: test if graph is acyclic (with TSort)
 * =====================================================================================
 */
c_bit** transitiveReduction(graph* p_gr)
{
  /* ORIGINAL ACCESSIBILITY MATRIX WITH 0'S -- no loops */
  int i, j;
  c_bit** orig = (c_bit **) calloc(p_gr->size, sizeof(c_bit *));
  for (i = 0; i < p_gr->size; i++) {
    orig[i] = (c_bit *) calloc(p_gr->size, sizeof(c_bit));
    for (j = 0; j < p_gr->size; j++)
      orig[i][j].bit = 0;
  }

  arc* p_cur_arc;
  for (i = 0; i < p_gr->size; i++) {
    p_cur_arc = p_gr->nodes[i].successors;       /* neighbors of i-th node */

    while (NULL != p_cur_arc) {
      orig[i][p_cur_arc->index].bit = 1;
      p_cur_arc = p_cur_arc->next;
    }
  }

  /* MATRIX OF TRANSITIVE CLOSURE */
  c_bit** closure = transitiveClosure(p_gr);

  /* reset the diagonal to 0's */
  for ( i = 0; i < p_gr->size; i += 1 ) {
    closure[i][i].bit = 0;
  }

  /* MATRIX OF (R ° R+) COMPOSITION */
  c_bit** composition = (c_bit **) calloc(p_gr->size, sizeof(c_bit *));
  for ( i = 0; i < p_gr->size; i += 1 ) {
    composition[i] = (c_bit *) calloc(p_gr->size, sizeof(c_bit));
    for ( j = 0; j < p_gr->size; j += 1 )
      composition[i][j].bit = 0;
  }

  /* fill the relation composition */
  int k;
  for ( i = 0; i < p_gr->size; i += 1 ) {
    for ( j = 0; j < p_gr->size; j += 1 ) {
      if ( 1 == closure[i][j].bit ) {             /* then copy j-th row to composition */
        for ( k = 0; k < p_gr->size; k += 1 ) {
          composition[i][k].bit |= orig[j][k].bit;
        }
      }
    }
  }

  c_bit** res = (c_bit **) calloc(p_gr->size, sizeof(c_bit *));
  for ( i = 0; i < p_gr->size; i += 1 ) {
    res[i] = (c_bit *) calloc(p_gr->size, sizeof(c_bit));
    for ( j = 0; j < p_gr->size; j += 1 ) {
      res[i][j].bit = orig[i][j].bit & ~composition[i][j].bit;
    }
  }

  int tmp_seq[p_gr->size];
  if ( -1 == TSort(p_gr, tmp_seq) ) {
    printf("Warning! Given graph is NOT acyclic.\n");
    printf("The result of transitive reduction is ungranted!\n");
  }
  return res;
}		/* -----  end of function transitiveReduction  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  matrixConvert
 *  Description:  receives pointer to list of successor, gives out adjacancy matrix
 * =====================================================================================
 */
double ** matrixConvert(graph* p_gr) {
  /* INIT OF MATRIX WITH 0'S */
  int i, j;
  double** matrix = (double **) calloc(p_gr->size, sizeof(double *));
  for (i = 0; i < p_gr->size; i++) {
    matrix[i] = (double *) calloc(p_gr->size, sizeof(double));
    for (j = 0; j < p_gr->size; j++)
      matrix[i][j] = DBL_MAX;
  }

  /* FILL MATRIX WITH ACTUAL WEIGHTS */
  arc* p_cur_arc;
  for (i = 0; i < p_gr->size; i++) {
    p_cur_arc = p_gr->nodes[i].successors;              /* neighbors of i-th node */
    
    while (NULL != p_cur_arc) {
      matrix[i][p_cur_arc->index] = p_cur_arc->weight;
      p_cur_arc = p_cur_arc->next;
    }
  }

  return matrix;
}		/* -----  end of matrixConvert  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  listConvert
 *  Description:  receives 2-D array of (double) weights of arcs, gives out the
 *  list-of-successors representation
 * =====================================================================================
 */
graph* listConvert(double** matrix, int dim)
{
  graph* p_gr = (graph*) malloc(sizeof(graph));
  p_gr->size = dim;
  p_gr->order = 0;
  p_gr->nodes = (node*) calloc(dim, sizeof(node));

  int i, j;
  for ( i = 0; i < dim; i += 1 ) {
    p_gr->nodes[i].successors = NULL;
  }

  arc* p_new_arc;
  for ( i = 0; i < dim; i += 1 ) {
    for ( j = 0; j < dim; j += 1 ) {
      if (DBL_MAX != matrix[i][j] && !(i == j && 0 == matrix[i][j])) {
        p_new_arc = (arc*) malloc(sizeof(arc));
        p_new_arc->index = j;
        p_new_arc->weight = matrix[i][j];
        p_new_arc->next = p_gr->nodes[i].successors;

        p_gr->nodes[i].successors = p_new_arc;
        p_gr->order++;
      }
    }
  }

  return p_gr;
}		/* -----  end of function listConvert  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  listConvertBit
 *  Description:  receives 2-D array of (double) weights of arcs, gives out the
 *  list-of-successors representation
 * =====================================================================================
 */
graph* listConvertBit(c_bit** matrix, int dim)
{
  graph* p_gr = (graph*) malloc(sizeof(graph));
  p_gr->size = dim;
  p_gr->order = 0;
  p_gr->nodes = (node*) calloc(dim, sizeof(node));

  int i, j;
  for ( i = 0; i < dim; i += 1 ) {
    p_gr->nodes[i].successors = NULL;
  }

  arc* p_new_arc;
  for ( i = 0; i < dim; i += 1 ) {
    for ( j = 0; j < dim; j += 1 ) {
      if (1 == matrix[i][j].bit) {
        p_new_arc = (arc*) malloc(sizeof(arc));
        p_new_arc->index = j;
        p_new_arc->weight = 1;
        p_new_arc->next = p_gr->nodes[i].successors;

        p_gr->nodes[i].successors = p_new_arc;
        p_gr->order++;
      }
    }
  }

  return p_gr;
}
