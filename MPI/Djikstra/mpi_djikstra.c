/* File:     mpi_dijkstra.c

   Author: Lyndon Ong Yiu
 * Purpose:  Implement Dijkstra's algorithm for solving the single-source 
 *           shortest path problem:  find the length of the shortest path 
 *           between a specified vertex and all other vertices in a 
 *           directed graph, dividing the work among multiple processes
 *           for efficiency
 *
 * Input:    n, the number of vertices in the digraph
 *           mat, the adjacency matrix of the digraph
 * Output:   A list showing the cost of the shortest path
 *           from vertex 0 to every other vertex in the graph.
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * Compile:  mpicc -g -Wall -o mstra mpi_dijkstra.c
 * Run:      mpiexec -n x ./mstra
 *           
 *           Note: The run statement includes variable x
 *           to represent number of processors. x should 
 *           be an int which can evenly divide n, number
 *           of vertices such that n % x = 0
 * *
 *
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#define MAX_STRING 10000
#define INFINITY 1000000

int Read_n(int my_rank, MPI_Comm comm);
MPI_Datatype Build_blk_col_type(int n, int loc_n);
void Read_matrix(int loc_mat[], int n, int loc_n, 
MPI_Datatype blk_col_mpi_t, int my_rank, MPI_Comm comm);
void Print_dists(int dist[], int n, int my_rank, int loc_n, MPI_Comm comm);
void Print_paths(int pred[], int loc_n, int my_rank, int n, MPI_Comm comm);
int  Find_min_dist(int dist[], int known[], int n, int my_rank, MPI_Comm comm);
void Dijkstra(int mat[], int loc_dist[], int loc_pred[], int loc_n, int n, int my_rank, MPI_Comm comm);

int main(int argc, char* argv[]) {
   int *loc_mat, *loc_dist, *loc_pred;
   int n, loc_n, p, my_rank;
   MPI_Comm comm;
   MPI_Datatype blk_col_mpi_t;
   MPI_Init(&argc, &argv);
   comm = MPI_COMM_WORLD;
   MPI_Comm_size(comm, &p);
   MPI_Comm_rank(comm, &my_rank);

   n = Read_n(my_rank, comm);
   loc_n = n/p;
   blk_col_mpi_t = Build_blk_col_type(n, loc_n);
   loc_mat = malloc(n*loc_n*sizeof(int));
   loc_pred = malloc(loc_n*sizeof(int));
   loc_dist = malloc(loc_n*sizeof(int));

   Read_matrix(loc_mat, n, loc_n, blk_col_mpi_t, my_rank, comm);
   Dijkstra(loc_mat, loc_dist, loc_pred, loc_n, n, my_rank, comm);
   Print_dists(loc_dist, n, my_rank, loc_n, comm);
   Print_paths(loc_pred, loc_n, my_rank, n, comm);

   free(loc_dist);
   free(loc_pred);
   free(loc_mat);
   MPI_Type_free(&blk_col_mpi_t);
   MPI_Finalize();

   return 0;
}  /* main */



/*-------------------------------------------------------------------
 * Function:    Dijkstra
 * Purpose:     Apply Dijkstra's algorithm to each process's loc_mat
 * In args:     loc_mat:  processor's local vertex distance values
 *              loc_dist: local distances to vertex v
 *              loc_pred:  local predecessors to vertex v
 *              loc_n:    number of vertices divded amongst processors
 *              n:        number of verteces
 *              my_rank:  processor rank
 *              MPI_Comm comm: MPI Communicator  
 *                
 * Out args:    dist:  loc_dist[v] = distance 0 to v.
 *              pred:  loc_pred[v] = predecessor of v on a 
 *                  shortest path 0->v.
 */
void Dijkstra(int loc_mat[], int loc_dist[], int loc_pred[], int loc_n, int n, int my_rank, MPI_Comm comm) {
   int i, v, best, loc_u, g_u, new_dist;
   int* loc_known;
   loc_known = malloc(loc_n*sizeof(int));
   for (v = 0; v < loc_n; v++) {
      loc_dist[v] = loc_mat[v];
      loc_pred[v] = 0;
      loc_known[v] = 0;
   }
   if (my_rank == 0){
      loc_known[0] = 1;
   }
   for (i = 1; i < n; i++) {
      g_u = Find_min_dist(loc_dist, loc_known, loc_n, my_rank, comm);
      if (g_u == -1){
         break;
      }
      loc_u = g_u % loc_n;
      if (g_u / loc_n == my_rank){
        
         loc_known[loc_u] = 1;
         best = loc_dist[loc_u];
      }

      MPI_Bcast(&best, 1, MPI_INT, g_u/loc_n, comm);
      for (v = 0; v < loc_n; v++){
         if (!loc_known[v]) {
            new_dist = best + loc_mat[g_u*loc_n + v];
            if (new_dist < loc_dist[v]) {
               loc_dist[v] = new_dist;
               loc_pred[v] = g_u;
            }
         }
      }
   } /* for i */
   free(loc_known);
}  /* Dijkstra */



/*-------------------------------------------------------------------
 * Function:    Find_min_dist
 * Purpose:     Find the vertex u with minimum distance to 0
 *              (dist[u]) among the vertices whose distance 
 *              to 0 is not known.
 * In args:     loc_dist:  loc_dist[v] = current estimate of distance
 *                 0->v
 *              loc_known:  whether the minimum distance 0-> is
 *                 known
 *              loc_n:  the total number of vertices per processor
 *              my_rank:  processor rank
 *              MPI_Comm comm: MPI Communicator  
 *
 * Ret val:     The vertex u whose distance to 0, loc_dist[u]
 *              is a minimum among vertices whose distance
 *              to 0 is not known.
 */
int Find_min_dist(int loc_dist[], int loc_known[], int loc_n, int my_rank, MPI_Comm comm) {
   int v, u, g_u, best_so_far = INFINITY;
   int g_min[2];
   int my_min[2];
   u = -1;
   for (v = 1; v < loc_n; v++){
      if (!loc_known[v]){
         if (loc_dist[v] < best_so_far) {
            u = v;
            best_so_far = loc_dist[v];
         }
      }
   }
   g_u = u + my_rank * loc_n;
   my_min[0] = best_so_far;
   my_min[1] = g_u;
   MPI_Allreduce(my_min, g_min, 1, MPI_2INT, MPI_MINLOC, comm);

   return g_min[1];
}  /* Find_min_dist */


/*-------------------------------------------------------------------
 * Function:    Print_dists
 * Purpose:     Print the length of the shortest path from 0 to each
 *              vertex
 * In args:     loc_dist: local distances to vertex v
 *              n:        number of verteces
 *              my_rank:  processor rank
 *              loc_n:    number of vertices divded amongst processors
 *              MPI_Comm comm: MPI Communicator  
 *                
 */
void Print_dists(int loc_dist[], int n, int my_rank, int loc_n, MPI_Comm comm) {
   int v;
   int* dist = NULL;
   if (my_rank ==0) dist = malloc(n *sizeof(int));
   MPI_Gather(loc_dist, loc_n, MPI_INT, dist, loc_n, MPI_INT, 0, comm);
   if (my_rank==0){
      printf("  v    dist 0->v\n");
      printf("----   ---------\n");          
      for (v = 1; v < n; v++){
         printf("%3d       %4d\n", v, dist[v]);
      }
      printf("\n");
      free(dist);
   }
} /* Print_dists */  


/*-------------------------------------------------------------------
 * Function:    Print_paths
 * Purpose:     Print the shortest path from 0 to each vertex
 * In args:     loc_pred: local predecessors to vertex 
 *              n:        number of vertices
 *              my_rank:  processor rank
 *              loc_n:    number of vertices divded amongst processors
 *              MPI_Comm comm: MPI Communicator  
 *                
 */
void Print_paths(int loc_pred[], int loc_n, int my_rank, int n, MPI_Comm comm) {
   int v, w, count, i;
   int* path = NULL;
   int* pred = NULL;
   if (my_rank == 0){
      path = malloc(loc_n*sizeof(int));
      pred = malloc(loc_n*sizeof(int));
   }
   MPI_Gather(loc_pred, loc_n, MPI_INT, pred, loc_n, MPI_INT, 0, comm);
   if (my_rank == 0){
      printf("  v     Path 0->v\n");
      printf("----    ---------\n");
      for (v = 1; v < n; v++) {
         printf("%3d:    ", v);
         count = 0;
         w = v;
         while (w != 0) {
            path[count] = w;
            count++;
            w = pred[w];
         }
         printf("0 ");
         for (i = count-1; i >= 0; i--){ 
            printf("%d ", path[i]);
         }
         printf("\n");
      }
      free(path);
      free(pred);
   }
   
}  /* Print_paths */

/*---------------------------------------------------------------------
 * Function:  Read_n
 * Purpose:   Read in the number of rows in the matrix on process 0
 *            and broadcast this value to the other processes
 * In args:   my_rank:  the calling process' rank
 *            comm:  Communicator containing all calling processes
 * Ret val:   n:  the number of rows in the matrix
 */
int Read_n(int my_rank, MPI_Comm comm) {
   int n;

   if (my_rank == 0){
      printf("Enter number of vertices\n");
      scanf("%d", &n);
   }
   MPI_Bcast(&n, 1, MPI_INT, 0, comm);
   return n;
}  /* Read_n */

/*---------------------------------------------------------------------
 * Function:  Read_matrix
 * Purpose:   Read in an nxn matrix of ints on process 0, and
 *            distribute it among the processes so that each
 *            process gets a block column with n rows and n/p
 *            columns
 * In args:   n:  the number of rows in the matrix and the submatrices
 *            loc_n = n/p:  the number of columns in the submatrices
 *            blk_col_mpi_t:  the MPI_Datatype used on process 0
 *            my_rank:  the caller's rank in comm
 *            comm:  Communicator consisting of all the processes
 * Out arg:   loc_mat:  the calling process' submatrix (needs to be 
 *               allocated by the caller)
 */
void Read_matrix(int loc_mat[], int n, int loc_n, 
      MPI_Datatype blk_col_mpi_t, int my_rank, MPI_Comm comm) {
   int* mat = NULL, i, j;

   if (my_rank == 0) {
      mat = malloc(n*n*sizeof(int));
      printf("Enter %d integers into the matrix\n", n*n);
      for (i = 0; i < n; i++){
         for (j = 0; j < n; j++){
            scanf("%d", &mat[i*n + j]);
         }
      }
   }

   MPI_Scatter(mat, 1, blk_col_mpi_t,
           loc_mat, n*loc_n, MPI_INT, 0, comm);

   if (my_rank == 0){
       free(mat);
   }
}  /* Read_matrix */


/*---------------------------------------------------------------------
 * Function:  Build_blk_col_type
 * Purpose:   Build an MPI_Datatype that represents a block column of
 *            a matrix
 * In args:   n:  number of rows in the matrix and the block column
 *            loc_n = n/p:  number cols in the block column
 * Ret val:   blk_col_mpi_t:  MPI_Datatype that represents a block
 *            column
 */
MPI_Datatype Build_blk_col_type(int n, int loc_n) {
   MPI_Aint lb, extent;
   MPI_Datatype block_mpi_t;
   MPI_Datatype first_bc_mpi_t;
   MPI_Datatype blk_col_mpi_t;

   MPI_Type_contiguous(loc_n, MPI_INT, &block_mpi_t);
   MPI_Type_get_extent(block_mpi_t, &lb, &extent);

   MPI_Type_vector(n, loc_n, n, MPI_INT, &first_bc_mpi_t);
   MPI_Type_create_resized(first_bc_mpi_t, lb, extent,
         &blk_col_mpi_t);
   MPI_Type_commit(&blk_col_mpi_t);

   MPI_Type_free(&block_mpi_t);
   MPI_Type_free(&first_bc_mpi_t);

   return blk_col_mpi_t;
}  /* Build_blk_col_type */
