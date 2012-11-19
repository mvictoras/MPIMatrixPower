/**
* Copyright (c) <2012>, <Victor Mateevitsi> www.vmateevitsi.com
* All rights reserved.
* 
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
* 3. All advertising materials mentioning features or use of this software
*    must display the following acknowledgement:
*    This product includes software developed by the <organization>.
* 4. Neither the name of the <organization> nor the
*    names of its contributors may be used to endorse or promote products
*    derived from this software without specific prior written permission.

* THIS SOFTWARE IS PROVIDED BY <COPYRIGHT HOLDER> ''AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
*	ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
*	(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
*	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*/

#include <getopt.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

int n, k, *local_coords, n_local;
float one_prob, minus_one_prob;
char s_local_coords[255];
int computerStats = 0;
int *counts;
int *displs;
MPI_Datatype resizedtype;
MPI_Comm comm_hor, comm_vert;

// Timing
double gen_time, proc_time, comm_time, total_time;
int hasReceivedAkk = -1;
typedef enum { 
  serial,
  strassen,
  cannon,
  dns
} TYPE;
TYPE type;

int parse_arguments(int argc, char **argv);
float *generate_array(int num_procs, char *proc_name, int local_rank, float one_prob, float minus_one_prob);
void create_2dmesh_topology(MPI_Comm *comm_new, int *local_rank, int *num_procs);
void create_ring_topology(MPI_Comm *comm_new, int *local_rank, int *num_procs);
void create_hypercube_topology(MPI_Comm *comm_new, int *local_rank, int *num_procs);

float *generate_serial(MPI_Comm *comm_new, int local_rank, int num_procs, char *proc_name, int *elem_per_node);
float *generate_parallel(MPI_Comm *comm_new, int local_rank, int num_procs, char *proc_name, int *elem_per_node);
float *distribute_hypercube(MPI_Comm *comm_new, int local_rank, int num_procs, char *proc_name, int *elem_per_node);

void serial_deter();

void determinant(float *A, int _n);
void one_d_partitioning(MPI_Comm *comm_new, float *A, int local_rank, int num_procs);
void two_d_partitioning(MPI_Comm *comm_new, float *A, int local_rank, int num_procs);

float *do_cannon(MPI_Comm *comm_new, float *A, float *B, int local_rank, int num_procs);
void do_strassen(float *A, float *B, float *C, int _n); 
void do_dns(MPI_Comm *comm_new, float *A, float *B, float *C, int local_rank, int num_procs);


void add(float *A, float *B, float *C, int _n);  
void sub(float *A, float *B, float *C, int _n);

void matrix_power(float *A, float *C, int _n, int _k);
void matrix_mult_serial(float *A, float *B, float *C, int n);
void block_mat_mult(float *A, float *B, float *C, int q);

void process_row_and_column(float *A, float *left_row, float *top_row, int startingRow, int k, int endingRow, int startingColumn, int endingColumn, int numRows, int numColumns, int local_k);

void printMatrix(float *A, int nElem);
void send_to(MPI_Comm *comm, int direction, float *A, int size, int row, int column, int n);
void receive_from_left(MPI_Comm *comm, int direction, float *A, int size, int row, int column, int n, int k);

int main(int argc, char **argv) {
	double t_start, t_end;
  int *ptr_gen_array, elem_per_node;
  float *C, *local_array;  
  int i, num_procs, local_rank, name_len;
  double start, end, dt;

	char proc_name[MPI_MAX_PROCESSOR_NAME];
	MPI_Comm comm_new;
  gen_time = 0.0; proc_time = 0.0; comm_time = 0.0; total_time = 0.0;
 
  // Parse the arguments
  if( parse_arguments(argc, argv) ) return 1;
  // Initialize MPI
  MPI_Init(&argc, &argv); 
	MPI_Get_processor_name(proc_name, &name_len);
  
  // Initially create the topology
  if( type == serial || type == strassen ) {
    create_ring_topology(&comm_new, &local_rank, &num_procs);
  } else if( type == cannon ) {
    create_2dmesh_topology(&comm_new, &local_rank, &num_procs);
  } else if( type == dns ) {
    create_hypercube_topology(&comm_new, &local_rank, &num_procs);
  }

  t_start = MPI_Wtime();
  if( type == serial ) {
    local_array = generate_serial(&comm_new, local_rank, num_procs, proc_name, &elem_per_node);
    //printMatrix(local_array, n);
    C = (float*) malloc(sizeof(float) * n * n);
    //bzero(C, sizeof(float) * n * n);
    start = MPI_Wtime();
    matrix_power(local_array, C, n, k);
    end = MPI_Wtime();
    dt = end - start;
    proc_time += dt;


    //printf("======\n");
    //printMatrix(C, n);
    start = MPI_Wtime();
    determinant(C, n);
    end = MPI_Wtime();
    dt = end - start;
    proc_time += dt;


    free(C);
  } else if( type == cannon ) {
    local_array = generate_parallel(&comm_new, local_rank, num_procs, proc_name, &elem_per_node);
    float *B, *C, *result;
    B = (float *)malloc(sizeof(float) * n_local * n_local);
    result = (float *)malloc(sizeof(float) * n * n);
    memcpy(B, local_array, sizeof(float) * n_local * n_local);
    
    for( i = 0; i < k - 1; ++i) {
      C = do_cannon(&comm_new, local_array, B, local_rank, num_procs);
      memcpy(local_array, C, sizeof(float) * n_local * n_local);
    }

    MPI_Gatherv(C, n_local * n_local, MPI_FLOAT, result, counts, displs, resizedtype, 0, comm_new); 
   
    if( local_rank == 0 ) {
      determinant(result, n);
    }
    free(C);
    free(B);
    free(result);
    MPI_Type_free(&resizedtype);
  
  } else if( type == strassen ) {
    float *B;
    local_array = generate_serial(&comm_new, local_rank, num_procs, proc_name, &elem_per_node);
    C = (float*)malloc(sizeof(float) * n * n);
    B = (float*)malloc(sizeof(float) * n * n);
    memcpy(B, local_array, sizeof(float) * n * n);

    for( i = 0; i < k - 1; ++i ) {
      do_strassen(local_array, B, C, n);
      memcpy(local_array, C, sizeof(float) * n * n);
    }

    //printMatrix(B, n);
    //printMatrix(C, n);
    
    start = MPI_Wtime();
    determinant(C, n);
    end = MPI_Wtime();
    dt = end - start;
    proc_time += dt;


    free(B);
    free(C);
  } else if( type == dns ) {
    ////local_array = distribute_hypercube(&comm_new, local_rank, num_procs, proc_name, &elem_per_node);
  }

  /*
  if( type == ring ) {
    one_d_partitioning(&comm_new, local_array, local_rank, num_procs);
  } else {
    two_d_partitioning(&comm_new, local_array, local_rank, num_procs);
  }
*/
  t_end = MPI_Wtime();
  total_time = t_end - t_start;

  if( computerStats ) {
    printf("%d\tg\t%s\t%d\t%f\n", n, s_local_coords, num_procs, gen_time);
    printf("%d\tp\t%s\t%d\t%f\n", n, s_local_coords, num_procs, proc_time);
    printf("%d\tc\t%s\t%d\t%f\n", n, s_local_coords, num_procs, comm_time);
    printf("%d\tt\t%s\t%d\t%f\n", n, s_local_coords, num_procs, total_time);
  }

  free(local_array);
	free(counts);
  free(displs);
  
  if( type == dns ) {
    //MPI_Comm_free(&comm_hor);
    //MPI_Comm_free(&comm_vert);
  }

  //MPI_Comm_free(&comm_new);
	MPI_Finalize(); // Exit MPI
	return 0;
}

int parse_arguments(int argc, char **argv) {
	int c;
  int option_index = 0;
  static struct option long_options[] =
  {
    {"input_method1", required_argument,  0, 'q'},
    {"input_method2", required_argument,  0, 'w'},
    {0, 0, 0, 0}
  };

  char *result = NULL;
  char delims[] = "m";
	while( (c = getopt_long (argc, argv, "n:t:q:w:k:c", long_options, &option_index)) != -1 ) {
		switch(c) {
      case 'q':
        //printf("optarg:%s\n", optarg);
        result = strtok(optarg, delims);
        one_prob = atof(result);
        result = strtok(NULL, delims);
        minus_one_prob = atof(result);
        break;
      case 'w':
        break;
			case 'n':
				n = atoi(optarg);
				break;
      case 'k':
        k = atoi(optarg);
        break;
      case 'c':
        computerStats = 1;
        break;
      case 't':
        if( strcmp(optarg, "serial" ) == 0 ) type = serial;
        else if( strcmp(optarg, "cannon" ) == 0 ) type = cannon;
        else if( strcmp(optarg, "strassen" ) == 0 ) type = strassen;
        else if( strcmp(optarg, "dns" ) == 0 ) type = dns;
        else {
          fprintf( stderr, "Option -%c %s in incorrect. Allowed values are: serial, strassen, cannon, dns\n", optopt, optarg);
          return 1;
        }
        break;
      case '?':
				if( optopt == 'n' )
					fprintf (stderr, "Option -%c requires an argument.\n", optopt);
				else if (isprint (optopt))
					fprintf (stderr, "Unknown option `-%c'.\n", optopt);
				else
					fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
				return 1;
			default:
				fprintf(stderr, "Usage: %s -n <number of numbers> \n", argv[0]);
				fprintf(stderr, "\tExample: %s -n 1000\n", argv[0]);
				return 1;
		}
	}
	return 0;
}


void create_ring_topology(MPI_Comm *comm_new, int *local_rank, int *num_procs) {
	int dims[1], periods[1];
 
  MPI_Comm_size(MPI_COMM_WORLD, num_procs);

  dims[0] = *num_procs;
  periods[0] = 1;
  local_coords = (int *) malloc(sizeof(int) * 1);
  // Create the topology
	MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 0, comm_new);
	MPI_Comm_rank(*comm_new, local_rank);
	MPI_Comm_size(*comm_new, num_procs);
	MPI_Cart_coords(*comm_new, *local_rank, 1, local_coords);
  sprintf(s_local_coords, "[%d]", local_coords[0]);
}

void create_2dmesh_topology(MPI_Comm *comm_new, int *local_rank, int *num_procs) {
  int *dims, i, *periods, nodes_per_dim;
  
  MPI_Comm_size(MPI_COMM_WORLD, num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, local_rank);
 
  int dimension = 2;
  nodes_per_dim = (int) sqrt( (double) *num_procs );
  local_coords = (int *) malloc(sizeof(int) * dimension);
  dims = (int *) malloc(sizeof(int) * dimension);
  periods = (int *) malloc(sizeof(int) * dimension);
  for( i = 0; i < dimension; i++ ) {
    dims[i] = nodes_per_dim;
    periods[i] = 1;
  }

  MPI_Cart_create(MPI_COMM_WORLD, dimension, dims, periods, 0, comm_new);
	MPI_Comm_size(*comm_new, num_procs);
  MPI_Cart_coords(*comm_new, *local_rank, dimension, local_coords);
  sprintf(s_local_coords, "[%d][%d]", local_coords[0], local_coords[1]);
}

void create_hypercube_topology(MPI_Comm *comm_new, int *local_rank, int *num_procs) {
  int dims[3], i, periods[3];
  int keep_dims[3];
  MPI_Comm_size(MPI_COMM_WORLD, num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, local_rank);
  
  //int dimension = (int) (log(*num_procs) / log(2));
  int dimension = 3;
  local_coords = (int *) malloc(sizeof(int) * dimension);
  for( i = 0; i < dimension; i++ ) {
    dims[i] = (int) cbrt(*num_procs);
    periods[i] = 1;
  }

  MPI_Cart_create(MPI_COMM_WORLD, dimension, dims, periods, 0, comm_new);
	MPI_Comm_size(*comm_new, num_procs);
  MPI_Cart_coords(*comm_new, *local_rank, dimension, local_coords);
  s_local_coords[0] = '\0';
  for( i = 0; i < dimension; i++ ) {
    char number[10];
    sprintf(number, "[%d]", local_coords[i]);
    strcat(s_local_coords, number);
  }

  keep_dims[0] = 1; // i
  keep_dims[1] = 1; // j
  keep_dims[2] = 0; // k
  printf("%s: Bang!\n", s_local_coords);
  // Create a 2-D sub-topology for each of the k-ith dimension
  MPI_Cart_sub(*comm_new, keep_dims, &comm_hor);
  printf("%s: Big Bang!\n", s_local_coords);

  //keep_dims[0] = 1;
  //keep_dims[1] = 0;
  //keep_dims[2] = 1;
  //MPI_Cart_sub(*comm_new, keep_dims, &comm_vert);
}

float *generate_array(int num_procs, char *proc_name, int local_rank, 
                      float one_prob, float minus_one_prob) {
	unsigned int iseed = (unsigned int)time(NULL);
  int i,j; 
  float *gen_array;
  double start, end, dt;
 
  srand (iseed);
	gen_array = (float *)malloc(sizeof(float) * (n * n));
  
  start = MPI_Wtime();
  float a = 0;
  int c = 0;
  int d = 0;
  //int k = 0;
	for( i = 0; i < n; i++ ) {
    for( j = 0; j < n; j++ ) {
      double number = (double)rand() / RAND_MAX;
      //printf("%f\n", number);
      if( number <= one_prob ) { a = 1; c++; }
      else if( number - one_prob <= minus_one_prob ) { a = -1; d++; }
      else a = 0;
      //a = k;
      //k++;
		  gen_array[i * n + j] = a;
    }
  }
  end = MPI_Wtime();
  dt = end - start;
  gen_time = dt;
  int c_n = n * n * one_prob;
  int d_n = n * n * minus_one_prob;

  //printMatrix(gen_array, n);
  if( !computerStats )
    printf("(%s(%d/%d)%s: %d random numbers generated in %1.8fs one:%d/%d minus one:%d/%d\n", proc_name, local_rank, num_procs, s_local_coords, n * n, dt, c, c_n, d, d_n);
  
  return gen_array;
}


float *generate_serial(MPI_Comm *comm_new, int local_rank, int num_procs, 
                 char *proc_name, int *elem_per_node) {
  float *local_array;
  double start, end, dt;
  int i, j;

  if( local_rank == 0 ) {
    local_array = generate_array(num_procs, proc_name, local_rank, one_prob, minus_one_prob);
  } else local_array = (float *) malloc(sizeof(float) * n * n);
 
  if( !computerStats )
    printf("(%s(%d/%d)%s: It took %1.8fs to receive the sub-array\n", proc_name, local_rank, num_procs, s_local_coords, dt);

  return local_array;
}


float *generate_parallel(MPI_Comm *comm_new, int local_rank, int num_procs, 
                 char *proc_name, int *elem_per_node) {
  float *local_array;
  float *tmp_array;
  double start, end, dt;
  int i, j;
  int sq_num_procs = sqrt(num_procs);
  n_local = n / sq_num_procs;
  MPI_Status status;

  local_array = (float *) malloc(sizeof(float) * n_local * n_local);
  counts = (int *)malloc(sizeof(int) * num_procs);
  displs = (int *)malloc(sizeof(int) * num_procs);
  
  if( local_rank == 0 ) {
    tmp_array = generate_array(num_procs, proc_name, local_rank, one_prob, minus_one_prob);
    //printMatrix(tmp_array, n);
  }
  for( i = 0; i <num_procs; i++ ) {
    counts[i] = 1;
  }
  
  for( i = 0; i < sq_num_procs; i++ ) {
    for(j = 0; j < sq_num_procs; j++ ) {
      displs[i * sq_num_procs + j] = j + i * n; 
    }
  }

  MPI_Datatype type;
  int sizes[2]    = {n,n};  /* size of global array */
  int subsizes[2] = {n_local,n_local};  /* size of sub-region */
  int starts[2]   = {0,0}; 
  
  /* as before */
  MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_FLOAT, &type);  
  /* change the extent of the type */
  MPI_Type_create_resized(type, 0, n_local * sizeof(float), &resizedtype);
  MPI_Type_commit(&resizedtype);
  
  MPI_Scatterv(tmp_array, counts, displs, resizedtype, local_array, n_local * n_local, MPI_FLOAT, 0, *comm_new); 

  //printf("(%s(%d/%d)%s: %d\n", proc_name, local_rank, num_procs, s_local_coords, displs[local_rank]);
  //printMatrix(local_array, n_local);
  
  //MPI_Type_free(&resizedtype);
  MPI_Type_free(&type);

  if( local_rank == 0 ) {
    free(tmp_array);
  }
  if( !computerStats )
    printf("(%s(%d/%d)%s: It took %1.8fs to receive the sub-array\n", proc_name, local_rank, num_procs, s_local_coords, dt);

  return local_array;
}

float *distribute_hypercube(MPI_Comm *comm_new, int local_rank, int num_procs, 
                 char *proc_name, int *elem_per_node) {
  
  float *local_array;
  float *tmp_array;
  double start, end, dt;
  int i, j;
  int sq_num_procs = sqrt(num_procs);
  n_local = n / sq_num_procs;
  MPI_Status status;

  local_array = (float *) malloc(sizeof(float) * n_local * n_local);
  counts = (int *)malloc(sizeof(int) * num_procs);
  displs = (int *)malloc(sizeof(int) * num_procs);
  
  if( local_rank == 0 ) {
    printf("Generating array\n");
    tmp_array = generate_array(num_procs, proc_name, local_rank, one_prob, minus_one_prob);
    //printMatrix(tmp_array, n);
  }
  for( i = 0; i <num_procs; i++ ) {
    counts[i] = 1;
  }
  
  for( i = 0; i < sq_num_procs; i++ ) {
    for(j = 0; j < sq_num_procs; j++ ) {
      displs[i * sq_num_procs + j] = j + i * n; 
    }
  }

  MPI_Datatype type;
  int sizes[2]    = {n,n};  // size of global array 
  int subsizes[2] = {n_local,n_local};  // size of sub-region
  int starts[2]   = {0,0}; 
  
  // as before
  MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_FLOAT, &type);  
  // change the extent of the type
  MPI_Type_create_resized(type, 0, n_local * sizeof(float), &resizedtype);
  MPI_Type_commit(&resizedtype);
 
  //if( local_coords[2] == 0 )
    //MPI_Scatterv(tmp_array, counts, displs, resizedtype, local_array, n_local * n_local, MPI_FLOAT, 0, comm_hor); 

  //printf("(%s(%d/%d)%s: %d\n", proc_name, local_rank, num_procs, s_local_coords, displs[local_rank]);
  //printMatrix(local_array, n_local);
  
  //MPI_Type_free(&resizedtype);
  MPI_Type_free(&type);

  if( local_rank == 0 ) {
    free(tmp_array);
  }
  if( !computerStats )
    printf("(%s(%d/%d)%s: It took %1.8fs to receive the sub-array\n", proc_name, local_rank, num_procs, s_local_coords, dt);

  return local_array;
}


void matrix_power(float *A, float *C, int _n, int _k) {
  int i;
  float *tmp_res;
  tmp_res = (float*)malloc(sizeof(float) * _n * _n);
  memcpy(C, A, sizeof(float) * _n * _n);

  for( i = 0; i < _k - 1; i++ ) {
    matrix_mult_serial(C, A, tmp_res, _n);
    memcpy(C, tmp_res, sizeof(float) * _n * _n);
  }

  free(tmp_res);
}

void matrix_mult_serial(float *A, float *B, float *C, int _n) {
  bzero(C, sizeof(float) * _n * _n);
  block_mat_mult(A, B, C, _n);
}

void block_mat_mult(float *A, float *B, float *C, int _q) {
  int i, j, k;
  //bzero(C, sizeof(float) * _q * _q);
  for( i = 0; i < _q; i++ ) {
    for( j = 0; j < _q; j++ ) { 
      for( k = 0; k < _q; k++ ) {
        C[i * _q + j] += A[i * _q + k] * B[k * _q + j];
      }
    }
  }
}

void determinant(float *A, int _n) {
  int k, i, j;
  long double determinant;
  double start, end, dt;


  for( k = 0;  k < _n; k++ ) {
    for( j = k + 1; j < _n; j++ ) {
      A[k * _n + j] = A[k * _n + j] / A[k * _n + k];
    }
    for( i = k + 1; i < _n; i++ ) {
      for( j = k + 1; j < _n; j++ ) {
        A[i * _n + j] -= A[i * _n + k] * A[k * _n + j];
      }
      A[i * _n + k] = 0;
    }
  }
  
  determinant = 1.0f;
  for( i = 0; i < n; i++ ) {
    //printf("%f " , A[i * n + i]);
    determinant = determinant * A[i * n + i];
  }
  //printf("\nDet is: %Lf\n", determinant);
}


void one_d_partitioning(MPI_Comm *comm_new, float *A, int local_rank, int num_procs) {
  MPI_Status status;
  int k, i, j, startingRow, endingRow, numRows;
  long double determinant;
  double start, end, dt;

  numRows = n / num_procs;
  startingRow = local_rank * numRows;
  endingRow = startingRow + numRows;

  start = MPI_Wtime();
  if( local_rank != 0 ) {
    MPI_Recv(A, n * n, MPI_FLOAT, local_rank - 1, 0, *comm_new, &status);  
  }
  end = MPI_Wtime();
  dt = end - start;
  comm_time += dt;
 
  start = MPI_Wtime();
  for( k = startingRow;  k < endingRow; k++ ) {
    for( j = k + 1; j < n; j++ ) {
      A[k * n + j] = A[k * n + j] / A[k * n + k];
    }
    for( i = k + 1; i < n; i++ ) {
      for( j = k + 1; j < n; j++ ) {
        A[i * n + j] -= A[i * n + k] * A[k * n + j];
      }
      A[i * n + k] = 0;
    }
  }
  end = MPI_Wtime();
  dt = end - start;
  proc_time += dt;

  start = MPI_Wtime();
  if( local_rank != num_procs - 1 )
    MPI_Send(A, n * n, MPI_FLOAT, local_rank + 1, 0, *comm_new);
  
  end = MPI_Wtime();
  dt = end - start;
  comm_time += dt;

  determinant = 1.0f;
  if( !computerStats && local_rank == num_procs - 1) {
   
    for( i = 0; i < n; i++ ) {
      //printf("%f " , A[i * n + i]);
      determinant = determinant * A[i * n + i];
    }
    //printf("\nDet is: %Lf\n", determinant);
  }
  end = MPI_Wtime();
  dt = end - start;
  proc_time += dt;
}

void two_d_partitioning(MPI_Comm *comm_new, float *A, int local_rank, int num_procs) {
  MPI_Status status;
  int k, i, j, startingRow, endingRow, numRows, startingColumn, endingColumn, numColumns;
  int n_startingRow, n_startingColumn, n_local_coords[2];
  //long double determinant;
  double start, end, dt;
  int p = (int) sqrt(num_procs);
  int dis, left_rank, right_rank, up_rank, down_rank;
  MPI_Request req;
  
  numRows = n / p;
  numColumns = numRows;

  startingRow = local_coords[1] * numRows;
  endingRow = startingRow + numRows;

  startingColumn = local_coords[0] * numRows;
  endingColumn = startingColumn + numColumns;

  
  start = MPI_Wtime();
  
  for( k = 0; k < n; k++ ) {
    float Akk[1];
    int local_k = k % numRows;
    
    // Send A(k,k) to the right
    start = MPI_Wtime();
    if( k >= startingColumn && k < endingColumn && k >= startingRow && k < endingRow ) {
      send_to(comm_new, 0, A, 1, local_k, local_k, numRows);
      Akk[0] = A[local_k * numRows + local_k];
    } else if( k < startingColumn && k >= startingRow && k < endingRow ) {
      receive_from_left(comm_new, 0, Akk, 1, 0, 0, numRows, k);
    }
    end = MPI_Wtime();
    dt = end - start;
    comm_time += dt;

    // Now calculate the row
    start = MPI_Wtime();
    if( k >= startingColumn && k < endingColumn && k >= startingRow && k < endingRow ) {
      for( j = local_k + 1; j < numColumns; j++ ) {
        A[local_k * numRows + j] /= Akk[0];
      }
    } else if( k >= startingRow && k < endingRow && k < startingColumn ) {
      for( j = 0; j < numColumns; j++ ) {
        A[local_k * numRows + j] /= Akk[0];
      }
    }
    end = MPI_Wtime();
    dt = end - start;
    proc_time += dt;
printf("%d: ASDAASD\n", local_rank);
    // Now calculate the box
    int m, bOutside = 1; 
    float top_row[numRows]; 

    start = MPI_Wtime();
    // k is West of this Partition
    if( k >= startingRow && k < endingRow & k < startingColumn ) {
      send_to(comm_new, 1, A, numColumns, local_k, 0, numRows);
      
      for( m = 0; m < numColumns; m++ ) {
        top_row[m] = A[local_k * numRows + m];
      }
      bOutside = -1;
    } 
    // k is in this BOX
    else if( k >= startingRow && k < endingRow && k >= startingColumn && k < endingColumn ) {
      int size = numColumns - (local_k + 1);
      if( size != 0 ) {
        send_to(comm_new, 1, A, size, local_k, local_k + 1, numRows);

        for( m = 0; m < size; m++ ) {
          top_row[m] = A[local_k * numRows + local_k + 1 + m];
        }
        bOutside = -1;
      }
    } // k is NW of this box 
    else if( k < startingRow && k < startingColumn ) {
      int sender_row = k / numRows;
      int sender_column = k / numColumns;
      int sender_rank = local_coords[0] * sqrt(num_procs) + sender_row;
      
      MPI_Recv(top_row, numColumns, MPI_FLOAT, sender_rank, 0, *comm_new, &status);
      
      bOutside = -1;
    } 
    // k is N of this box
    else if( k < startingRow && k >= startingColumn && k < endingColumn ) {
      int sender_row = k / numRows;
      int sender_column = k / numColumns;
      int sender_rank = sender_column * sqrt(num_procs) + sender_row;
      int size = numColumns - (local_k + 1);
      
      if( size != 0 ) { 
        //top_row = (float *)malloc(sizeof(float) * numberToReceive);
        //printf("%d Waiting to receive from:%d\n", local_rank, sender_rank);
        MPI_Recv(top_row, size, MPI_FLOAT, sender_rank, 0, *comm_new, &status);
        
        bOutside = -1;
      } 
    }
    float left_row[numRows];
    
    // k is N of this Box
    if( k >= startingColumn && k < endingColumn & k < startingRow ) {
      for(m = 0; m < numRows; m++ ) {
        left_row[m] = A[m * numColumns + local_k];
      }
      send_to(comm_new, 0, left_row, numRows, 0, 0, 0);
      
      bOutside = -1;
    } 
    // k is IN this box 
    else if( k >= startingRow && k < endingRow && k >= startingColumn && k < endingColumn ) {
      //int local_k = k % numRows;
      int size = numColumns - (local_k + 1);
      if( size != 0 ) {
        for(m = 0; m < size; m++ ) {
          left_row[m] = A[(local_k + 1) * numColumns + local_k];
        }
        send_to(comm_new, 0, left_row, size, 0, 0, 0);
        
        bOutside = -1;
      }
    } 
    // k is SW from this box
    else if( k < startingRow && k < startingColumn ) {
      int sender_row = k / numRows;
      int sender_column = k / numColumns;
      int sender_rank = sender_column * sqrt(num_procs) + local_coords[1];
      
      MPI_Recv(left_row, numColumns, MPI_FLOAT, sender_rank, 0, *comm_new, &status);
      
      bOutside = -1;
    } 
    // k is W of this box
    else if( k < startingColumn && k >= startingRow && k < endingRow ) {
      int sender_row = k / numRows;
      int sender_column = k / numColumns;
      int sender_rank = sender_column * sqrt(num_procs) + local_coords[1];
      int local_k = k % numRows;
      int numberToReceive = numColumns - (local_k + 1);
      
      if( numberToReceive != 0 ) { 
        MPI_Recv(left_row, numberToReceive, MPI_FLOAT, sender_rank, 0, *comm_new, &status);
        bOutside = -1;
      }
    }
    end = MPI_Wtime();
    dt = end - start;
    comm_time += dt;
   

    // Now process the box
    if( bOutside < 0 ) {
      start = MPI_Wtime();
      process_row_and_column(A, left_row, top_row, k, startingRow, endingRow, startingColumn, endingColumn, numRows, numColumns, local_k);
      end = MPI_Wtime();
      dt = end - start;
      proc_time += dt;
    }
  } // end for

  float determinant[1];
  float result[1];
  determinant[0] = 1;
  if( local_coords[0] == local_coords[1] ) {
    start = MPI_Wtime();
    for(i = 0; i < numRows; i++ ) {
      determinant[0] *= A[i * numRows + i];
    }
    end = MPI_Wtime();
    dt = end - start;
    proc_time += dt;
  }
  
  start = MPI_Wtime();
  MPI_Reduce(determinant, result, 1, MPI_FLOAT, MPI_PROD, 0, *comm_new);
  end = MPI_Wtime();
  dt = end - start;
  comm_time += dt;
  if( !computerStats && local_rank == 0 ) {
    printf("Determinant is %f\n", result[0]);
  }
}

void do_strassen(float *A, float *B, float *C, int _n) {
  int i, j;

  int half_n = _n / 2;
  float *a11, *a12, *a21, *a22, *b11, *b12, *b21, *b22;
  float *c11, *c12, *c21, *c22;
  float *aRes, *bRes;
  float *m1, *m2, *m3, *m4, *m5, *m6, *m7;

  if( _n == 1 ) {
    C[0] = A[0] * B[0];
    return;
  }

  // Allocate the memories
  a11 = (float *)malloc(sizeof(float) * half_n * half_n);
  a12 = (float *)malloc(sizeof(float) * half_n * half_n);
  a21 = (float *)malloc(sizeof(float) * half_n * half_n);
  a22 = (float *)malloc(sizeof(float) * half_n * half_n);

  b11 = (float *)malloc(sizeof(float) * half_n * half_n);
  b12 = (float *)malloc(sizeof(float) * half_n * half_n);
  b21 = (float *)malloc(sizeof(float) * half_n * half_n);
  b22 = (float *)malloc(sizeof(float) * half_n * half_n);

  aRes = (float *)malloc(sizeof(float) * half_n * half_n);
  bRes = (float *)malloc(sizeof(float) * half_n * half_n);

  m1 = (float *)malloc(sizeof(float) * half_n * half_n);
  m2 = (float *)malloc(sizeof(float) * half_n * half_n);
  m3 = (float *)malloc(sizeof(float) * half_n * half_n);
  m4 = (float *)malloc(sizeof(float) * half_n * half_n);
  m5 = (float *)malloc(sizeof(float) * half_n * half_n);
  m6 = (float *)malloc(sizeof(float) * half_n * half_n);
  m7 = (float *)malloc(sizeof(float) * half_n * half_n);

  for( i = 0; i < half_n; i++ ) {
    for( j = 0; j < half_n; j++ ) {
      a11[i * half_n + j] = A[i * _n + j];
      a12[i * half_n + j] = A[i * _n + j + half_n];
      a21[i * half_n + j] = A[(i + half_n) * _n + j];
      a22[i * half_n + j] = A[(i + half_n) * _n + j + half_n];
      
      b11[i * half_n + j] = B[i * _n + j];
      b12[i * half_n + j] = B[i * _n + j + half_n];
      b21[i * half_n + j] = B[(i + half_n) * _n + j];
      b22[i * half_n + j] = B[(i + half_n) * _n + j + half_n];
    }
  }

  add(a11, a22, aRes, half_n); 
  add(b11, b22, bRes, half_n);
  do_strassen(aRes, bRes, m1, half_n); // m1 = (a11+a22) * (b11+b22)
                   
  add(a21, a22, aRes, half_n); 
  do_strassen(aRes, b11, m2, half_n); // m2 = (a21+a22) * (b11)
                             
  sub(b12, b22, bRes, half_n);
  do_strassen(a11, bRes, m3, half_n); // m3 = (a11) * (b12 - b22)
                                                   
  sub(b21, b11, bRes, half_n);
  do_strassen(a22, bRes, m4, half_n); // m4 = (a22) * (b21 - b11)
                                                                 
  add(a11, a12, aRes, half_n);
  do_strassen(aRes, b22, m5, half_n); // m5 = (a11+a12) * (b22)   
                                                                                  
  sub(a21, a11, aRes, half_n);
  add(b11, b12, bRes, half_n); 
  do_strassen(aRes, bRes, m6, half_n); // m6 = (a21-a11) * (b11+b12)
                                                         
  sub(a12, a22, aRes, half_n);
  add(b21, b22, bRes, half_n);
  do_strassen(aRes, bRes, m7, half_n); // m7 = (a12-a22) * (b21+b22)

  free(a11);
  free(a12);
  free(a21);
  free(a22);

  free(b11);
  free(b12);
  free(b21);
  free(b22);

  c11 = (float *)malloc(sizeof(float) * half_n * half_n);
  c12 = (float *)malloc(sizeof(float) * half_n * half_n);
  c21 = (float *)malloc(sizeof(float) * half_n * half_n);
  c22 = (float *)malloc(sizeof(float) * half_n * half_n);

  // Calculating C matrices
  add(m3, m5, c12, half_n); // c12 = p3 + p5
  add(m2, m4, c21, half_n); // c21 = p2 + p4
  
  add(m1, m4, aRes, half_n); // p1 + p4
  add(aRes, m7, bRes, half_n); // p1 + p4 + p7
  sub(bRes, m5, c11, half_n); // c11 = p1 + p4 - p5 + p7
                                   
  add(m1, m3, aRes, half_n); // p1 + p3
  add(aRes, m6, bRes, half_n); // p1 + p3 + p6
  sub(bRes, m2, c22, half_n); // c22 = p1 + p3 - p2 + p6
  
  for( i = 0; i < half_n; ++i ) {
    for( j = 0 ; j < half_n; ++j ) {
      C[i * _n + j]                     = c11[i * half_n + j];
      C[i * _n + j + half_n]            = c12[i * half_n + j];
      C[(i + half_n) * _n + j]          = c21[i * half_n + j];
      C[(i + half_n) * _n + j + half_n] = c22[i * half_n + j];
    }
  }

  free(c11);
  free(c12);
  free(c21);
  free(c22);

  free(aRes);
  free(bRes);

  free(m1);
  free(m2);
  free(m3);
  free(m4);
  free(m5);
  free(m6);
  free(m7);
}

float *do_cannon(MPI_Comm *comm_new, float *A, float *B, int local_rank, int num_procs) {
  float *C;
  int i;
  MPI_Status status;
  int left_rank, right_rank, up_rank, down_rank;
  int distance = 1;
  double start, end, dt;

  //for( i = 0; i < num_procs; ++i ) {
  //  MPI_Barrier( *comm_new );
  //  if ( i == local_rank ) {
      //printf("B:%d:\n", local_rank);
      //printMatrix(A, n);
  //  }
 // }
 
  C = (float *) malloc(sizeof(float) * n_local * n_local);
  bzero(C, sizeof(float) * n_local * n_local);
  //printMatrix(A, n);
  // Initial alignement
  MPI_Cart_shift(*comm_new, 1, distance, &left_rank, &right_rank);
  MPI_Cart_shift(*comm_new, 0, distance, &up_rank, &down_rank);
 
  start = MPI_Wtime();
  for( i = 0; i < local_coords[0]; i++ ) {
    MPI_Sendrecv_replace(A, n_local * n_local, MPI_FLOAT, left_rank, 0, right_rank, 0, *comm_new, &status);
  }
  
  for( i = 0; i < local_coords[1]; i++ ) {
    MPI_Sendrecv_replace(B, n_local * n_local, MPI_FLOAT, up_rank, 0, down_rank, 0, *comm_new, &status);
  }
  end = MPI_Wtime();
  dt = end - start;
  comm_time += dt;

  
  start = MPI_Wtime();
  block_mat_mult(A, B, C, n_local);
  end = MPI_Wtime();
  dt = end - start;
  proc_time += dt;


  /*for( i = 0; i < num_procs; ++i ) {
    MPI_Barrier( *comm_new );
    if ( i == 0 && local_rank == 0 ) {
      printf("A:%d:\n", local_rank);
      printMatrix(A, n_local);
      printMatrix(B, n_local);
      printMatrix(C, n_local);
    }
  }*/
  for( i = 0; i < sqrt(num_procs) - 1; ++i ) {
    //matrix_power(A, C, n_local, k);
    //MPI_Cart_shift(*comm_new, 1, distance, &left_rank, &right_rank);
    start = MPI_Wtime();
    MPI_Sendrecv_replace(A, n_local * n_local, MPI_FLOAT, left_rank, 0, right_rank, 0, *comm_new, &status);
    
    //MPI_Cart_shift(*comm_new, 0, distance, &up_rank, &down_rank);
    MPI_Sendrecv_replace(B, n_local * n_local, MPI_FLOAT, up_rank, 0, down_rank, 0, *comm_new, &status);
    end = MPI_Wtime();
    dt = end - start;
    comm_time += dt;

//printf("rank:%d A:%f B:%f\n", local_rank, A[0], B[0]);
    start = MPI_Wtime();
    block_mat_mult(A, B, C, n_local);
    end = MPI_Wtime();
    dt = end - start;
    proc_time += dt;


  //printf("rank:%d up:%d down:%d left:%d right%d\n", local_rank, up_rank, down_rank, left_rank, right_rank);
  }

  
  //for( i = 0; i < num_procs; ++i ) {
  //  MPI_Barrier( *comm_new );
  //  if ( i == local_rank ) {
  //    printf("A:%d:\n", local_rank);
      //printMatrix(A, n_local);
      //printMatrix(B, n_local);
  //    printMatrix(C, n_local);
  //  }
  //}

  return C;
}

void  do_dns(MPI_Comm *comm_new, float *A, float *B, float *C, int local_rank, int num_procs) {
}

void process_row_and_column(float *A, float *left_row, float *top_row, int k,
                            int startingRow, int endingRow, int startingColumn, int endingColumn, 
                            int numRows, int numColumns, int local_k) {
  int i, j;
  int index_row = 0;
  int index_column = 0;
  int starting_x = 0;
  int starting_y = 0;
  if( k >= startingColumn && k < endingColumn ) {
    starting_y = local_k + 1;
  }
  if( k >= startingRow && k < endingRow ) {
    starting_x = local_k + 1;
  }

  char str[255];
  for( i = starting_x; i < numRows; i++ ) {
    index_column = 0;
    for( j = starting_y; j < numColumns; j++ ) {
      A[i * numRows + j] -= left_row[index_row] * top_row[index_column];    
      index_column++;;
    }
    index_row++;
  }
}

void printMatrix(float *A, int nElem) {
  int i, j;
  for( i = 0; i < nElem; i++ ) {
    for( j = 0; j < nElem; j++ ) {
      printf("%f ", A[i * nElem + j]);
    }
    printf("\n");
  }
  printf("\n");
}

/*
 * int direction (0 horizontal, 1 vertical)
 * int distance  
 */
void send_to(MPI_Comm *comm, int direction, float *A, int size, int row, int column, int n) {
  int prev_rank, next_rank;
  int distance = 1;
  MPI_Cart_shift(*comm, direction, distance, &prev_rank, &next_rank);
  while(next_rank >= 0) {
    MPI_Send(A + row * n + column, size, MPI_FLOAT, next_rank, 0, *comm);
    MPI_Cart_shift(*comm, direction, ++distance, &prev_rank, &next_rank);
  }
} 
    
void receive_from_left(MPI_Comm *comm, int direction, float *A, int size, int row, int column, int n, int k) {
  int prev_rank, next_rank, coords[2], startingRow, startingColumn;
  int distance = 1;
  MPI_Status status;

  MPI_Cart_shift(*comm, direction, distance, &prev_rank, &next_rank);
  while(prev_rank >= 0) {
    MPI_Cart_coords(*comm, prev_rank, 2, coords);
    
    startingRow     = coords[1] * n; 
    startingColumn  = coords[0] * n;
    
    if( k >= startingColumn && k < startingColumn + n ) {
      MPI_Recv(A, size, MPI_FLOAT, prev_rank, 0, *comm, &status);
    }
    MPI_Cart_shift(*comm, 0, ++distance, &prev_rank, &next_rank);
  }
}

void add(float *A, float *B, float *C, int _n) {  
  int i, j;
  
  for (i = 0; i < _n; ++i) {
    for (j = 0; j < _n; ++j) {
      C[i * _n + j] = A[i * _n + j] + B[i * _n + j];
    }
  }

}

void sub(float *A, float *B, float *C, int _n) {  
  int i, j;
  
  for (i = 0; i < _n; ++i) {
    for (j = 0; j < _n; ++j) {
      C[i * _n + j] = A[i * _n + j] - B[i * _n + j];
    }
  }

}

