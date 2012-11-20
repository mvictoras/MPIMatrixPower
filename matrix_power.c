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
int assoc = 0;
int input_numbers[4];
int input_method2 = 0;

MPI_Datatype resizedtype;
MPI_Comm comm_hor, comm_vert, comm_vert_col, comm_hor_row, comm_hor_col;

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

float *distribute_serial(MPI_Comm *comm_new, int local_rank, int num_procs, char *proc_name, int *elem_per_node);
float *distribute_mesh(MPI_Comm *comm_new, int local_rank, int num_procs, char *proc_name, int *elem_per_node);
void distribute_hypercube(MPI_Comm *comm_new, int local_rank, int num_procs, char *proc_name, int *elem_per_node, float *A, float *B);

void determinant(float *A, int _n);

float *do_cannon(MPI_Comm *comm_new, float *A, float *B, int local_rank, int num_procs);
void do_strassen(float *A, float *B, float *C, int _n); 
void do_dns(MPI_Comm *comm_new, float *A, float *B, float *C, int local_rank, int num_procs);

void dns_align(float *A, float *B);

void add(float *A, float *B, float *C, int _n);  
void sub(float *A, float *B, float *C, int _n);

void matrix_power(float *A, float *C, int _n, int _k);
void matrix_mult_serial(float *A, float *B, float *C, int n);
void block_mat_mult(float *A, float *B, float *C, int q);


void printMatrix(float *A, int nElem);

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
    local_array = distribute_serial(&comm_new, local_rank, num_procs, proc_name, &elem_per_node);
    C = (float*) malloc(sizeof(float) * n * n);
    start = MPI_Wtime();
    matrix_power(local_array, C, n, k);
    end = MPI_Wtime();
    dt = end - start;
    proc_time += dt;


    start = MPI_Wtime();
    determinant(C, n);
    end = MPI_Wtime();
    dt = end - start;
    proc_time += dt;


    free(C);
  } else if( type == cannon ) {
    local_array = distribute_mesh(&comm_new, local_rank, num_procs, proc_name, &elem_per_node);
    float *B, *C, *result;
    B = (float *)malloc(sizeof(float) * n_local * n_local);
    result = (float *)malloc(sizeof(float) * n * n);
    memcpy(B, local_array, sizeof(float) * n_local * n_local);
    
    int k_end = (assoc) ? k / 2 : k;
    
    for( i = 2; i <= k_end; ++i) {
      C = do_cannon(&comm_new, local_array, B, local_rank, num_procs);
      memcpy(local_array, C, sizeof(float) * n_local * n_local);
    }

    if( assoc ) {
      float *Bcopy = (float *) malloc(sizeof(float) * n_local * n_local);
      memcpy(Bcopy, B, sizeof(float) * n_local * n_local);
      
      memcpy(B, local_array, sizeof(float) * n_local * n_local);
      C = do_cannon(&comm_new, local_array, B, local_rank, num_procs);
      
      if( k % 2 != 0 ) {
        memcpy(local_array, C, sizeof(float) * n_local * n_local);
        memcpy(B, Bcopy, sizeof(float) * n_local * n_local);
        C = do_cannon(&comm_new, local_array, B, local_rank, num_procs);
      }
      free(Bcopy);
    }

    start = MPI_Wtime();
    MPI_Gatherv(C, n_local * n_local, MPI_FLOAT, result, counts, displs, resizedtype, 0, comm_new); 
    end = MPI_Wtime();
    dt = end - start;
    proc_time += dt;

    if( local_rank == 0 ) {
      start = MPI_Wtime();
      determinant(result, n);
      end = MPI_Wtime();
      dt = end - start;
      proc_time += dt;
    }
    free(C);
    free(B);
    free(result);
    MPI_Type_free(&resizedtype);
  
  } else if( type == strassen ) {
    float *B;
    local_array = distribute_serial(&comm_new, local_rank, num_procs, proc_name, &elem_per_node);
    C = (float*)malloc(sizeof(float) * n * n);
    B = (float*)malloc(sizeof(float) * n * n);
    memcpy(B, local_array, sizeof(float) * n * n);

    int k_end = (assoc) ? k / 2 : k;
 
    for( i = 2; i <= k_end; ++i ) {
      do_strassen(local_array, B, C, n);
      memcpy(local_array, C, sizeof(float) * n * n);
    }
   
    if( assoc ) {
      float *Bcopy = (float *) malloc(sizeof(float) * n * n);
      memcpy(Bcopy, B, sizeof(float) * n * n);
      
      memcpy(B, local_array, sizeof(float) * n * n);
      do_strassen(local_array, B, C, n);
      
      if( k % 2 != 0 ) {
        memcpy(local_array, C, sizeof(float) * n * n);
        memcpy(B, Bcopy, sizeof(float) * n * n);
        do_strassen(local_array, B, C, n);
      }
      free(Bcopy);
    }
    
    start = MPI_Wtime();
    determinant(C, n);
    end = MPI_Wtime();
    dt = end - start;
    proc_time += dt;


    free(B);
    free(C);
  } else if( type == dns ) {
    float *A, *B, *C, *result;
    int cbrt_num_procs = cbrt(num_procs);
    n_local = n / cbrt_num_procs;
    
    A = (float *) malloc(sizeof(float) * n_local * n_local);
    B = (float *) malloc(sizeof(float) * n_local * n_local);
    C = (float *) malloc(sizeof(float) * n_local * n_local);
    result = (float *) malloc(sizeof(float) * n * n);

    distribute_hypercube(&comm_new, local_rank, num_procs, proc_name, &elem_per_node, A, B);

    // Align B
    MPI_Bcast(B, n_local * n_local, MPI_FLOAT, 0, comm_vert_col);
    MPI_Bcast(B, n_local * n_local, MPI_FLOAT, local_coords[2], comm_hor_col);

    int k_end = (assoc) ? k / 2 : k;
    
    for( i = 2; i <= k_end; ++i) {
      // Align A
      MPI_Bcast(A, n_local * n_local, MPI_FLOAT, 0, comm_vert_col);
      MPI_Bcast(A, n_local * n_local, MPI_FLOAT, local_coords[2], comm_hor_row);
      do_dns(&comm_new, A, B, C, local_rank, num_procs);
      memcpy(A, C, sizeof(float) * n_local * n_local);
    }

    if( assoc ) {
      float *Bcopy = (float *) malloc(sizeof(float) * n_local * n_local);
      memcpy(Bcopy, B, sizeof(float) * n_local * n_local);
      
      memcpy(B, A, sizeof(float) * n_local * n_local);
      // Align B 
      MPI_Bcast(B, n_local * n_local, MPI_FLOAT, 0, comm_vert_col);
      MPI_Bcast(B, n_local * n_local, MPI_FLOAT, local_coords[2], comm_hor_col);
      
      // Align A
      MPI_Bcast(A, n_local * n_local, MPI_FLOAT, 0, comm_vert_col);
      MPI_Bcast(A, n_local * n_local, MPI_FLOAT, local_coords[2], comm_hor_row);
      do_dns(&comm_new, A, B, C, local_rank, num_procs);
      if( k % 2 != 0 ) {
        memcpy(A, C, sizeof(float) * n_local * n_local);
        memcpy(B, Bcopy, sizeof(float) * n_local * n_local);

        MPI_Bcast(A, n_local * n_local, MPI_FLOAT, 0, comm_vert_col);
        MPI_Bcast(A, n_local * n_local, MPI_FLOAT, local_coords[2], comm_hor_row);
        do_dns(&comm_new, A, B, C, local_rank, num_procs);
      }
      free(Bcopy);
    }
    
    MPI_Gatherv(C, n_local * n_local, MPI_FLOAT, result, counts, displs, resizedtype, 0, comm_hor); 

    start = MPI_Wtime();
    determinant(result, n);
    end = MPI_Wtime();
    dt = end - start;
    proc_time += dt;

    free(A);
    free(B);
    free(C);
    free(result);
  }

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
    MPI_Comm_free(&comm_vert_col);
    MPI_Comm_free(&comm_hor_row);
    MPI_Comm_free(&comm_hor_col);
    MPI_Comm_free(&comm_hor);
    MPI_Comm_free(&comm_vert);
  }

  MPI_Comm_free(&comm_new);
	MPI_Finalize(); // Exit MPI
	return 0;
}

int parse_arguments(int argc, char **argv) {
	int i, c;
  int option_index = 0;
  static struct option long_options[] =
  {
    {"input_method1", required_argument,  0, 'q'},
    {"input_method2", required_argument,  0, 'w'},
    {0, 0, 0, 0}
  };

  char *result = NULL;
  char delims[] = "m";
	while( (c = getopt_long (argc, argv, "n:t:q:w:k:ca", long_options, &option_index)) != -1 ) {
		switch(c) {
      case 'q':
        result = strtok(optarg, delims);
        one_prob = atof(result);
        result = strtok(NULL, delims);
        minus_one_prob = atof(result);
        break;
      case 'w':
        input_method2 = 1;
        result = strtok(optarg, delims);
        input_numbers[0] = atoi(result);
        for( i = 1; i < 4; i++ ) {
          result = strtok(NULL, delims);
          input_numbers[i] = atoi(result);
        }
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
      case 'a':
        assoc = 1;
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

  // Create a 2-D sub-topology for each of the k-ith dimension
  keep_dims[0] = 1; // i
  keep_dims[1] = 1; // j
  keep_dims[2] = 0; // k
  MPI_Cart_sub(*comm_new, keep_dims, &comm_hor);

  // Create one along the j dimension
  keep_dims[0] = 1;
  keep_dims[1] = 0;
  keep_dims[2] = 1;
  MPI_Cart_sub(*comm_new, keep_dims, &comm_vert);

  keep_dims[0] = 0;
  keep_dims[1] = 0;
  keep_dims[2] = 1;
  MPI_Cart_sub(*comm_new, keep_dims, &comm_vert_col);

  keep_dims[0] = 0;
  keep_dims[1] = 1;
  keep_dims[2] = 0;
  MPI_Cart_sub(*comm_new, keep_dims, &comm_hor_row);

  keep_dims[0] = 1;
  keep_dims[1] = 0;
  keep_dims[2] = 0;
  MPI_Cart_sub(*comm_new, keep_dims, &comm_hor_col);

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
  int input_i = 0;
	
  for( i = 0; i < n; i++ ) {
    for( j = 0; j < n; j++ ) {
      if( !input_method2 ) {
        double number = (double)rand() / RAND_MAX;
        if( number <= one_prob ) { a = 1; c++; }
        else if( number - one_prob <= minus_one_prob ) { a = -1; d++; }
        else a = 0;
      } else {
        a = input_numbers[input_i];
        input_i = (input_i == 3) ? 0 : input_i + 1; 
      }
		  gen_array[i * n + j] = a;
    }
  }
  end = MPI_Wtime();
  dt = end - start;
  gen_time = dt;
  int c_n = n * n * one_prob;
  int d_n = n * n * minus_one_prob;

  if( !computerStats )
    printf("(%s(%d/%d)%s: %d random numbers generated in %1.8fs one:%d/%d minus one:%d/%d\n", proc_name, local_rank, num_procs, s_local_coords, n * n, dt, c, c_n, d, d_n);
  
  return gen_array;
}


float *distribute_serial(MPI_Comm *comm_new, int local_rank, int num_procs, 
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


float *distribute_mesh(MPI_Comm *comm_new, int local_rank, int num_procs, 
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
 
  start = MPI_Wtime();
  MPI_Scatterv(tmp_array, counts, displs, resizedtype, local_array, n_local * n_local, MPI_FLOAT, 0, *comm_new); 
  end = MPI_Wtime();
  dt = end - start;
  comm_time += dt;

  MPI_Type_free(&type);

  if( local_rank == 0 ) {
    free(tmp_array);
  }
  if( !computerStats )
    printf("(%s(%d/%d)%s: It took %1.8fs to receive the sub-array\n", proc_name, local_rank, num_procs, s_local_coords, dt);

  return local_array;
}

void distribute_hypercube(MPI_Comm *comm_new, int local_rank, int num_procs, 
                 char *proc_name, int *elem_per_node, float *A, float *B) {
  
  float *tmp_array;
  double start, end, dt;
  int i, j;
  int cbrt_num_procs = cbrt(num_procs);
  MPI_Status status;

  counts = (int *)malloc(sizeof(int) * num_procs);
  displs = (int *)malloc(sizeof(int) * num_procs);
  
  if( local_rank == 0 ) {
    tmp_array = generate_array(num_procs, proc_name, local_rank, one_prob, minus_one_prob);
  }
  for( i = 0; i <num_procs; i++ ) {
    counts[i] = 1;
  }
  
  for( i = 0; i < cbrt_num_procs; i++ ) {
    for(j = 0; j < cbrt_num_procs; j++ ) {
      displs[i * cbrt_num_procs + j] = j + i * n; 
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

  start = MPI_Wtime();
  if( local_coords[2] == 0 ) {
    MPI_Scatterv(tmp_array, counts, displs, resizedtype, A, n_local * n_local, MPI_FLOAT, 0, comm_hor); 
    MPI_Scatterv(tmp_array, counts, displs, resizedtype, B, n_local * n_local, MPI_FLOAT, 0, comm_hor); 
  }
  end = MPI_Wtime();
  dt = end - start;
  comm_time += dt;

  
  MPI_Type_free(&type);

  if( local_rank == 0 ) {
    free(tmp_array);
  }
}

void dns_align(float *A, float *B) {
  // Align A
  MPI_Bcast(A, n_local * n_local, MPI_FLOAT, 0, comm_vert_col);
  MPI_Bcast(A, n_local * n_local, MPI_FLOAT, local_coords[2], comm_hor_row);
  
  // Align B
  MPI_Bcast(B, n_local * n_local, MPI_FLOAT, 0, comm_vert_col);
  MPI_Bcast(B, n_local * n_local, MPI_FLOAT, local_coords[2], comm_hor_col);
}

void matrix_power(float *A, float *C, int _n, int _k) {
  int i;
  float *B;
  B = (float*)malloc(sizeof(float) * _n * _n);
  memcpy(B, A, sizeof(float) * _n * _n);

  int k_end = (assoc) ? _k / 2 : _k;
  for( i = 2; i <= k_end; ++i ) {
    matrix_mult_serial(A, B, C, _n);
    memcpy(A, C, sizeof(float) * _n * _n);
  }
  
  if( assoc ) {
    float *Bcopy = (float *) malloc(sizeof(float) * _n * _n);
    memcpy(Bcopy, B, sizeof(float) * _n * _n);
      
    memcpy(B, A, sizeof(float) * _n * _n);
    matrix_mult_serial(A, B, C, _n);
      
    if( k % 2 != 0 ) {
      memcpy(A, C, sizeof(float) * _n * _n);
      memcpy(B, Bcopy, sizeof(float) * _n * _n);
      matrix_mult_serial(A, B, C, _n);
    }
    free(Bcopy);
  }

  free(B);
}

void matrix_mult_serial(float *A, float *B, float *C, int _n) {
  bzero(C, sizeof(float) * _n * _n);
  block_mat_mult(A, B, C, _n);
}

void block_mat_mult(float *A, float *B, float *C, int _q) {
  int i, j, k;
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
    determinant = determinant * A[i * n + i];
  }
}

void do_strassen(float *A, float *B, float *C, int _n) {
  int i, j;
  double start, end, dt;
  start = MPI_Wtime();

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
  end = MPI_Wtime();
  dt = end - start;
  proc_time += dt;
}

float *do_cannon(MPI_Comm *comm_new, float *A, float *B, int local_rank, int num_procs) {
  float *C;
  int i;
  MPI_Status status;
  int left_rank, right_rank, up_rank, down_rank;
  int distance = 1;
  double start, end, dt;

  C = (float *) malloc(sizeof(float) * n_local * n_local);
  bzero(C, sizeof(float) * n_local * n_local);
  
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

  for( i = 0; i < sqrt(num_procs) - 1; ++i ) {
    start = MPI_Wtime();
    MPI_Sendrecv_replace(A, n_local * n_local, MPI_FLOAT, left_rank, 0, right_rank, 0, *comm_new, &status);
    
    MPI_Sendrecv_replace(B, n_local * n_local, MPI_FLOAT, up_rank, 0, down_rank, 0, *comm_new, &status);
    end = MPI_Wtime();
    dt = end - start;
    comm_time += dt;

    start = MPI_Wtime();
    block_mat_mult(A, B, C, n_local);
    end = MPI_Wtime();
    dt = end - start;
    proc_time += dt;
  }

  return C;
}

void  do_dns(MPI_Comm *comm_new, float *A, float *B, float *C, int local_rank, int num_procs) {
  int i;
  float *AB;
  double start, end, dt;
  AB = (float *) malloc(sizeof(float) * n_local * n_local);
  bzero(AB, sizeof(float) * n_local * n_local);

  start = MPI_Wtime();
  block_mat_mult(A, B, AB, n_local);
  end = MPI_Wtime();
  dt = end - start;
  proc_time += dt;

  start = MPI_Wtime();
  MPI_Reduce(AB, C, n_local * n_local, MPI_FLOAT, MPI_SUM, 0, comm_vert_col);
  end = MPI_Wtime();
  dt = end - start;
  comm_time += dt;

  free(AB);

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

