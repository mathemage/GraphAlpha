/*
 * main.c - GraphAlpha v1.3
 * =============================================================
 * Demonstration program of GraphAlpha (implementation of graph algorithms)
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include "GraphAlpha.h"

void displayMenu() {
  printf("---------------------------------------------------------------\n");
  printf("Choose option (enter the part in brackets)\n");
  printf("(S)ym \t\t - symmetrization\n");
  printf("(D)ijkstra \t - shortest path between node & other nodes\n");
  printf("(F)loyd-Warshall - shortest path between every 2 nodes\n");
  printf("tar(J)an \t - strongly connected components\n");
  printf("(T)sort \t - topological ordering\n");
  printf("(C)losure \t - transitive closure\n");
  printf("(R)eduction \t - transitive reduction\n");
  printf("l(I)st-convert \t - from adjacency matrix to list of successors\n");
  printf("(M)atrix-convert - from list of successors to adjacency matrix\n");
  printf("Displa(Y) \t - display the graph \n");
  printf("c(L)ear \t - clear screen (only for in *nix like OS) \n");
  printf("(H)elp \t\t - display this menu:) \n");
  printf("(Q)uit \t\t - quit the program \n");
  printf("---------------------------------------------------------------\n");
}

void displayGraph(graph* p_gr) {
  int i;
  arc* cur_arc;
  for (i = 0; i < p_gr->size; i++) {
    if (NULL == p_gr->nodes[i].successors) {        /* nothing-2-C */
      continue;
    }

    printf("%d. node ->", i+1);

    cur_arc = p_gr->nodes[i].successors; 
    while (NULL != cur_arc) {
      printf(" %d(%lf)", cur_arc->index+1, cur_arc->weight);
      cur_arc = cur_arc->next;
    }

    printf("\n");
  }
}

int main() {
  printf("Welcome to GraphAlpha v1.3!\n");
  printf("This is a demonstration program of GraphAlpha library.\n");

  /* LOADING GRAPH */
  graph gr;

  printf("Enter count of nodes: ");
  scanf("%d", &gr.size);
  gr.nodes = calloc(gr.size, sizeof(node));
  int i;
  for (i = 0; i < gr.size; i++) {      /* temporarily NULL pointers for successors */
    gr.nodes[i].successors = NULL;
  }

  printf("Enter count of arcs: ");
  scanf("%d", &gr.order);

  printf("Enter (directed) arcs\n");
  unsigned int start, end;
  double a_weight;
  arc* p_tmp_arc;
  for (i = 0; i < gr.order; i++) {
    printf("%d. node\n", i+1);

    do {
      printf(" From: "); scanf("%d", &start);
      printf(" To: "); scanf("%d", &end);
      printf(" Weight: "); scanf("%lf", &a_weight);
    } while (start < 1 || start > gr.size || end < 1 || end > gr.size);

    p_tmp_arc = (arc *) malloc(sizeof(arc));
    p_tmp_arc->index = end-1;
    p_tmp_arc->weight = a_weight;
    p_tmp_arc->next = gr.nodes[start-1].successors;
    gr.nodes[start-1].successors = p_tmp_arc;
  }

  displayMenu();
  
  scanf("%*c");
  char cOption;                         /* chosen command */
  do {
    printf(">> ");
    scanf("%c%*c", &cOption);           /* read the first char, omit the rest of line */
    cOption = toupper(cOption);

    /* what option? */
    switch (cOption) {
      /* symmetrization of arcs */
      case 'S': {
                  double** matrix;
                  matrix = symMtrx(&gr);

                  printf("List of successor representation:\n");
                  displayGraph(listConvert(matrix, gr.size));

                  printf("Matrix representation:\n");
                  int i, j;
                  for (i = 0; i < gr.size; i++) {
                    for (j = 0; j < gr.size; j++)
                      if (DBL_MAX == matrix[i][j])
                        printf("----\t");
                      else
                        printf("%.2lf\t", matrix[i][j]);
                    free(matrix[i]);
                    printf("\n");
                  }
                  free(matrix);
                }
                break;

      /* Dijkstra's algorithm shortest path between given node & the rest of graph */
      case 'D': {
                  printf("Enter starting node of graph: ");
                  int v;
                  scanf("%d%*c", &v);

                  double* distances = (double *) calloc(gr.size, sizeof(double));
                  distances = Dijkstra(v-1, &gr);
                  printf("From %d.node to all nodes respectively:\n", v);
                  int i;
                  for ( i = 0; i < gr.size; i += 1 ) {
                    if (DBL_MAX == distances[i]) {
                      printf("---- ");
                    }
                    else {
                      printf("%lf ", distances[i]);
                    }
                  }
                  printf("\n");
                }
                break;

      /* Floyd-Warshall's algorithm for finding shortest paths between every 2 nodes */
      case 'F': {
                  double** matrix;
                  matrix = FloydWarshall(&gr);

                  printf("List of successor representation");
                  printf(" (i.e., transitive closure):\n");
                  displayGraph(listConvert(matrix, gr.size));

                  printf("Matrix representation:\n");
                  int i, j;
                  for (i = 0; i < gr.size; i++) {
                    for (j = 0; j < gr.size; j++)
                      if (DBL_MAX == matrix[i][j])
                        printf("----\t");
                      else
                        printf("%.2lf\t", matrix[i][j]);
                    free(matrix[i]);
                    printf("\n");
                  }
                  free(matrix);
                }
                break;

      /* Tarjan's algorithm for strongly connected components */
      case 'J': Tarjan(&gr);
                break;

      /* topological ordering */
      case 'T': {
                  int seq[gr.size];
                  if (0 == TSort(&gr, seq) ) {
                    printf("Topological ordering: ");
                    for ( i = 0; i < gr.size; i += 1 ) {
                      printf("%d ", seq[i]+1);
                    }
                    printf("\n");
                  }
                  else {
                   printf("NOT AN ACYCLIC GRAPH!");
                  }
                }
                break;

      /* transitive closure */
      case 'C': {
                  c_bit** matrix;
                  matrix = transitiveClosure(&gr);

                  printf("List of successor representation:\n");
                  displayGraph(listConvertBit(matrix, gr.size));

                  printf("Matrix representation:\n");
                  int i, j;
                  for (i = 0; i < gr.size; i++) {
                    for (j = 0; j < gr.size; j++)
                      printf( (matrix[i][j].bit) ? "1 " : "0 ");
                    free(matrix[i]);
                    printf("\n");
                  }
                  free(matrix);
                }
                break;

      /* transitive reduction */
      case 'R': {
                  c_bit** matrix;
                  matrix = transitiveReduction(&gr);

                  printf("List of successor representation:\n");
                  displayGraph(listConvertBit(matrix, gr.size));

                  printf("Matrix representation:\n");
                  int i, j;
                  for (i = 0; i < gr.size; i++) {
                    for (j = 0; j < gr.size; j++)
                      printf( (matrix[i][j].bit) ? "1 " : "0 ");
                    free(matrix[i]);
                    printf("\n");
                  }
                  free(matrix);
                }
                break;

      /* convert from adjacency matrix to list of successors */
      case 'I': {
                  printf("Converting the already matrix-converted graph to ");
                  printf("list-of-sucessors representation...\n");

                  double** matrix = matrixConvert(&gr);
                  displayGraph((graph*) listConvert(matrix, gr.size));

                  int i;
                  for ( i = 0; i < gr.size; i += 1 ) {
                    free(matrix[i]);
                  }
                  free(matrix);
                }
                break;

      /* convert from list of successors to adjacency matrix */
      case 'M': {
                  double** matrix;
                  matrix = matrixConvert(&gr);
                  int i, j;
                  for (i = 0; i < gr.size; i++) {
                    printf("\n");
                    for (j = 0; j < gr.size; j++) {
                      if (DBL_MAX == matrix[i][j])
                        printf("----\t");
                      else
                        printf("%.2lf\t", matrix[i][j]);
                    }
                    free(matrix[i]);
                  }
                  free(matrix);
                }
                break;

      /* display command menu */
      case 'Y': displayGraph(&gr);
                break;

      /* display command menu */
      case 'H': displayMenu();
                break;

      /* clear screen - works only on *nix systems */
      case 'L': system("clear");
                break;

      /* delete the tree & quit */
      case 'Q': printf("Bye bye...\n");
                break;

      default:  printf("Invalid command \"%c\"! Please try again...\n", cOption);
    }
  } while ('Q' != cOption);

  return 0;
}
